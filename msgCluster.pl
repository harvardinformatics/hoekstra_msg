#!/usr/bin/perl -w
use strict;
use lib qw(./msg .);
use Utils;

print "\nMSG\n";
my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
  localtime(time);
printf "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year + 1900, $mon + 1, $mday, $hour,
  $min, $sec;

### Make sure all required dependencies are installed
&Utils::test_dependencies();

### Default parameters
### All of these parameters are required
my %default_params = (
	barcodes                            => 'NULL',
	re_cutter                           => 'MseI',
	linker_system                       => 'Dros_SR_vII',
	reads                               => 'NULL',
	parent1                             => 'NULL',
	parent2                             => 'NULL',
	chroms                              => 'all',
	sexchroms                           => 'X',
	chroms2plot                         => 'all',
	deltapar1                           => '.01',
	deltapar2                           => '.01',
	recRate                             => '0',
	rfac                                => '0.00001',
	thinfac                             => '1',
	difffac                             => '.01',
	priors                              => '0,.5,.5',
	bwaindex1                           => 'bwtsw',
	bwaindex2                           => 'bwtsw',
	pnathresh                           => '0.03',
	cluster                             => '1',
	threads                             => '8',
	theta                               => '1',
	addl_qsub_option_for_exclusive_node => '',
	addl_qsub_option_for_pe             => '',
	custom_qsub_options_for_all_cmds    => '',
	bwa_alg                             => 'aln',
	bwa_threads                         => '1',
	use_stampy                          => '0',
	stampy_premap_w_bwa                 => '1',
	stampy_pseudo_threads               => '0',
	quality_trim_reads_thresh           => '0',
	quality_trim_reads_consec           => '30',
	indiv_stampy_substitution_rate      => '0.001',
	indiv_mapq_filter                   => '0',
	index_file                          => '',
	index_barcodes                      => '',
	email_host                          => '',
	notify_emails                       => '',
	debug                               => '0',
	gff_thresh_conf                     => '.95',
	new_parser                          => '0',
	new_parser_offset                   => '0',
	new_parser_filter_out_seq           => '',
	pepthresh                           => '',
	one_site_per_contig                 => '0',
);

my $params = Utils::parse_config( 'msg.cfg', \%default_params );
Utils::validate_config( $params, qw( barcodes reads parent1 parent2 ) );
my %params = %$params;

### check if all the desired chroms are found in both parental files
### report their lengths also
my %par1_reads = &Utils::readFasta( $params{'parent1'}, 1 );
my %par2_reads = &Utils::readFasta( $params{'parent2'}, 1 );
my @chroms;
if ( $params{'chroms'} eq 'all' ) {
	@chroms = keys %par1_reads;
}
else { @chroms = split( /,/, $params{'chroms'} ); }

my $numcontigs = length(@chroms);

open( OUT, '>msg.chrLengths' )
  || die "ERROR (msgCluster): Can't create msg.chrLengths: $!\n";
print OUT "chr,length\n";
foreach my $chr ( sort @chroms ) { print OUT "$chr,$par1_reads{$chr}\n"; }
close OUT;

### Mapping & Plotting
### qsub array: one for each line in the barcode file
my $num_barcodes = 0;
open( FILE, $params{'barcodes'} )
  || die "ERROR (msgCluster): Can't open $params{'barcodes'}: $!\n";
while (<FILE>) {
	chomp $_;
	if ( $_ =~ /^\S+\t.*$/ ) {
		$num_barcodes++;
	}
}
close FILE;

print "num barcodes is $num_barcodes!\n";

####################################################################################################
mkdir "msgOut.$$"   unless ( -d "msgOut.$$" );
mkdir "msgError.$$" unless ( -d "msgError.$$" );

### Run jobs!
open (BARCODE,$barcodes) || die "ERROR: Can't open $barcodes: $!\n";
foreach my $bc_line (<BARCODE>) {
    chomp $bc_line;

    # fastq_file = 'indiv' + ind[1] + '_' + ind[0]
    my @bc_bits = split(/\s+/,$bc_line);
    my $indiv = 'indiv' . $bc_bits[1] . '_' . $bc_bits[0];
    print "\t$indiv\n";
    
    my $cmd="Rscript $src/fit-hmm.R -d $outdir -i $indiv -s $sex -o $Routdir -p $deltapar1 -q $deltapar2 -a $recRate -r $rfac -c $chroms -x $sexchroms -y $chroms2plot -z $priors -t $theta -g $gff_thresh_conf -u $one_site_per_contig -j $pepthresh"
    
}


&Utils::system_call(
	"python msg/create_stats.py -i $params{'reads'} -b $params{'barcodes'}"
);
if ( $params{'pepthresh'} ne '' ) {
	&Utils::system_call(
        "python msg/hmmprob_to_est.py -d hmm_fit -t $params{'pepthresh'} -o hmm_fits_ests.csv"
	);
}
&Utils::system_call(
    "Rscript msg/summaryPlots.R -c $params{'chroms'} -p $params{'chroms2plot'} -d hmm_fit -t $params{'thinfac'} -f $params{'difffac'} -b $params{'barcodes'} -n $params{'pnathresh'} > msgRun3.$$.out 2> msgRun3.$$.err"
);
&Utils::system_call(
    "perl msg/summary_mismatch.pl $params{'barcodes'} 0"
);

#Run a simple validation
&Utils::system_call(
    "python msg/validate.py $params{'barcodes'} > msgRun.validate.$$.out 2> msgRun.validate.$$.err"
);

#Cleanup - move output files to folders, remove barcode related files
&Utils::system_call(
    "mv -f msgRun*.${$}.e** msgError.$$; mv -f msgRun*.${$}.pe** msgError.$$; mv -f msgRun*.${$}.o* msgOut.$$; mv -f msgRun*.${$}.po* msgOut.$$; mv -f *.trim.log msgOut.$$; rm -f temp.fq; rm -f $params{'barcodes'}.*"
);

#Notify users that MSG run has completed
if ( $params{'email_host'} && $params{'notify_emails'} ) {
	&Utils::system_call( "python msg/send_email.py -e $params{'email_host'}"
		  . " -t $params{'notify_emails'} -s 'MSG Run has completed'"
		  . " -b 'NOTE: Output and error messages are located in: msgOut.$$ and msgError.$$'"
	);
}

print
"\nNOTE: Output and error messages are located in: msgOut.$$ and msgError.$$ \n\n";
exit;
