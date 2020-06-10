use strict;
use warnings;
use File::Temp qw/tempfile tempdir/;
use FindBin qw($Bin);
use POSIX;
use Getopt::Long;
use File::Basename;
use Env qw(@PATH);
use threads;
use Thread::Semaphore;
use threads::shared;
use File::Which;
use Statistics::R;
use Bio::DB::HTS;
use Statistics::Descriptive;
use Cwd qw(abs_path getcwd);

my $path = $Bin;
my $line;
my $bam;
my $out;
my $tmp;
my $bar_th;
my $seg_th;
my $gap;
my $size;
my $is;
my $Nmerge;
my $bin;
my $low;
my $ex;
my $phase_dir;
my $black_list;
my $control_list;
my $id_num;
my $Smerge;
my $ncpu;
my $human;
my $blow;
my $bex;
my $mergemax;
my $qc1;
my $qc2;
my $qc3;
my $rlen;
my $mlen;
my $help;


unshift @PATH, "$path/bin";
my $Basename=basename($0);
my $USAGE = qq{
Name:
	$Basename
	version 2.1.2
Function:
	Detect the SVs from stLFR WGS data
Usage:
	$Basename -bam prefix.sorted.markdup.bam -o out.result.dir
Options:
	-bam <string>	original sorted and markduped bam file,if the index dose not exist, will be created.[necessary]
	-out <string>	output SV dir.[necessary](warning:if exists, the output dir will be cleaned first!!!!!)
	-ncpu <int>	thread number for running pipeline.[default 1]
	-bar_th <int> at least N read pairs in one barcode.[default 8]
	-seg_th <int> at least N read pairs in one segment.[default 4]
	-gap <int> define the gap size which should be considered as different segment.[default 20000]
	-size <int> output SV length.[default 20000]
	-is <int> proper IS size for read pair library, read pairs with too large IS will be abandoned.[default 300]
	-bin <int> bin size for cluster the segments.[default 2000]
	-merge1 <int> N continue outline bins could be considered as the same break point and will be merged into one evidence.[default 5]
	-merge2 <int> SVs nearby under N binsize will be considered as one event.[default 5]
	-mmax <int> the max SVs allowed in one event.[default 4]
	-low1 <float/int> single end barcode counts threshold, 0-1 float: higher than X percentage counts; 1> int: higher than X counts.[default 0.95]
	-low2 <float/int> end to end barcode counts threshold, 0-1 float: higher than X percentage counts; 1> int: higher than X counts.[default 0.9995]
	-ex1 <float> when low1 is a float of 0-1, exclude the bins which depth under ex1.[default 0.2]
	-ex2 <float> when low2 is a float of 0-1, exclude the bins which depth under ex2.[default 0.2]
	-phase <string> formatted phase result directory including phased barcode and region by chromosome.[default NULL]
	-bl <string> black list file(BED format).[default NULL]
	-cl <string> sorted control list file(BEDPE format).[default NULL](Be sure the chromosome and position are sorted in one line!!!)
	-sc <int> allow max sv counts for the same position in one direction.[default 4]
	-human <Y/N> for Homo sapiens,keep only [1234567890XYM] chromosome.[default N]
	-qc1 <float> valid read pair ratio for SV detection.[default 0.60]
	-qc2 <float> average read pair count for one barcode.[default 30]
	-qc3 <float> average segment end count for one bin.[default 8]
	-rlen <int> read length of one read.[default 100]
	-mlen <int> physical limit for the long DNA segment.[default 400000]
	-help Show this message.
};

$ncpu ||= 1;
$tmp ||= "/tmp";
$bar_th||=8;
$seg_th ||= 4;
$gap ||= 20000;
$size ||= 20000;
$is ||=300;
$bin ||=2000;
$Nmerge ||=5;
$low ||=0.95;
$blow ||= 0.9995;
$ex ||=0.2;
$bex ||=0.2;
$phase_dir ||="NULL";
$black_list||="NULL";
$control_list||="NULL";
$id_num||=4;
$Smerge||=5;
$human||="N";
$mergemax ||=4;
$qc1 ||=0.6;
$qc2 ||=30;
$qc3 ||=8;
$rlen ||=100;
$mlen ||=400000;

my $valid;
$valid = GetOptions(
"help"=>\$help,
"bam=s"=>\$bam, 
"out=s"=>\$out,
"ncpu=i"=>\$ncpu,
"bar_th=i"=>\$bar_th,
"seg_th=i"=>\$seg_th,
"gap=i"=>\$gap,
"size=i"=>\$size,
"is=i" =>\$is,
"bin=i" =>\$bin,
"merge1=i" =>\$Nmerge,
"merge2=i" =>\$Smerge,
"low1=f" =>\$low,
"low2=f" =>\$blow,
"ex1=f" =>\$ex,
"ex2=f" =>\$bex,
"mmax=i" =>\$mergemax,
"phase=s" => \$phase_dir,
"bl=s" => \$black_list,
"cl=s" => \$control_list,
"sc=i" => \$id_num,
"human=s" => \$human,
"qc1=f" => \$qc1,
"qc2=f" => \$qc2,
"qc3=f" => \$qc3,
"rlen=i" => \$rlen,
"mlen=i" => \$mlen
);

die "$USAGE" unless (($bam and $out) and !$help);

unless(-d $out){
	`mkdir -p $out`;
	unless(-d $out){
		die "Cannot creat out dir $out\n";
	}
}else{
	`rm -rf $out`;
	`mkdir -p $out`;
	my @files=glob "$out/*";
	if(@files >0){
		die "Check the permission of out dir $out\n";
	}
}

which 'samtools' or die "samtools was not found\n";
which 'tabix' or die "tabix was not found\n";

if(!(-e $bam)){die "Cannot find the file $bam\n";}

$bam=abs_path($bam);
$out=abs_path($out);
my $bamname=basename($bam);
if(!(-e "$bam.bai")){
	print STDERR "Cannot find the index file $bam.bai\n";
	print STDERR "Will creat the bam index now\n";
	$line="ln -s $bam $out/";
	if(executeSystemCall($line)){
		die "Failure link the $bam\n";
	}
	$line="samtools index $out/$bamname";
	if(executeSystemCall($line)){
		die "Failure index the $bam\n";
	}
	$bam="$out/$bamname";
}

my $start = time;

print STDERR "Creat the sorted barcoded file\n";
$line="barcode-sort $bam $ncpu $out $bamname $bar_th $seg_th $gap $is $human";
if(executeSystemCall($line)){
	die "Failure during creating the sorted barcoded file $out/$bamname.sbf\n";
}
####QC#########
open SAM,"samtools view -H $bam|" or die $!;
my $reflen=0;
while(<SAM>){
	chomp;
	my $line=$_;
	next unless $line=~ /^\@SQ/;
	$line=~ /SN:(\S+)\s+LN:(\d+)/;
	my $name=$1;
	my $len=$2;
	if($human eq "Y"){
		next unless $name=~ /^(?:chr)?[1234567890XYM]{1,2}$/;
	}
	$reflen+=$len;
}
close SAM;

if($reflen == 0){
	die "Wrong bam header data for the reference!\nPlease check the Bam file...\n";
}

open IN,"$out/$bamname.stat" or die $!;
my $tl=<IN>;
close IN;
chomp $tl;
my ($tread,$vread,$vbar)=split(/\t/,$tl);
my $q1=$vread/$tread;
my $q2=$vread/$vbar;
my $q3=$vbar*2/($reflen/$bin);
if($q1 < $qc1){
	print STDERR "Warning: valid read pair ratio is lower than $qc1, the result may be unreliable!\n";
}

if($q2 < $qc2){
	print STDERR "Warning: average read pair on one barcode is lower than $qc2, the result may be unreliable!\n";
}

if($q3 < $qc3){
	print STDERR "Warning: average segment end in one bin is lower than $qc3, the result may be unreliable!\n";
}

my $R = Statistics::R->new();
my $cmd=qq{
a = read.table("$out/$bamname.gap")
x<-a[,1]
x<-x[which(x>$rlen)]
x<-x[which(x<$mlen)]
len1<-quantile(x,0.65)
len2<-quantile(x,0.93)
len3<-quantile(x,0.98)
};
$R->run($cmd);
my $s1=$R->get('len1');
my $s2=$R->get('len2');
my $s3=$R->get('len3');
$s2=int($s2/$s1+0.5);
print STDERR "According to the read pair distance statistic:\n";
print STDERR "bin size should be set around $s1\n";
print STDERR "merge bin number should be set around $s2\n";
print STDERR "gap size should be set around $s3\n";
###########################

print STDERR "Step_1 single end cluster\n";
$line="sin-cluster $out/$bamname.sbf $out/$bamname.sin $bin";

if(executeSystemCall($line)){
	die "Failure during creating the one end file $out/$bamname.sin\n";
}

print STDERR "Step_2 link the cluster to cluster\n";

$line="link-id $out/$bamname.sin $out/$bamname.lnd $low $ex $blow $bex $gap $size $bin $Nmerge";

if(executeSystemCall($line)){
	die "Failure during creating the link file $out/$bamname.lnd\n";
}

print STDERR "Step_3 relink the cluster to the segment\n";

$line="link-sin-seg $out/$bamname.sbf $out/$bamname.sin $out/$bamname.lns";

if(executeSystemCall($line)){
	die "Failure during creating the link file $out/$bamname.lns\n";
}

##the segment length is short, pass by the ins detect module
#split id link module
#step4
print STDERR "Step_4 split the link into haploidy\n";
$line="split-link $out/$bamname.sin $out/$bamname.lns $out/$bamname.lnd $out/$bamname.sln";

if(executeSystemCall($line)){
	die "Failure during creating the split file $out/$bamname.sln\n";
}

#step5
print STDERR "Step_5 judge the link Quality\n";
$line="judge-link $out/$bamname.sin $out/$bamname.sln $bam $phase_dir $out/$bamname.judge $bin $ncpu";
if(executeSystemCall($line)){
	die "Failure during creating the judgement file $out/$bamname.judge\n";
}

#step6
print STDERR "Step_6 filter by the LFR segment\n";
$line="id_filter $out/$bamname.judge $out/$bamname.sln $out/$bamname.filter $id_num";
if(executeSystemCall($line)){
	die "Failure during creating the filter file $out/$bamname.filter\n";
}


#step7
print STDERR "Step_7 filter by the Regions\n";
my $extend_len=$Nmerge*$bin;
$line="region_filter $out/$bamname.filter $black_list $control_list $out/$bamname.region $extend_len";
if(executeSystemCall($line)){
	die "Failure during creating the filter file $out/$bamname.region\n";
}

#step8
print STDERR "Step_8 link and merge SVs for final output\n";
my $merge_size=$Smerge*$bin;
$line="id_link $out/$bamname.region $out/$bamname.sln $out/$bamname.final $merge_size $mergemax";
if(executeSystemCall($line)){
	die "Failure during output the final file $out/$bamname.final\n";
}

#step9
print STDERR "Step_9 plot heatmaps for final output\n";
my $ext=2*$gap;
$line="bat_plot $bam $out/$bamname.final 0 $bin $ext 0 1 $out/heatmap_plot $ncpu 0";
if(executeSystemCall($line)){
	die "Failure during plotting the heatmaps\n";
}

my $end=time-$start;
print STDERR "Done. time: ",sprintf("%.2f",$end/3600)," hours\n";



sub executeSystemCall {
  my ($command,$returnVal) = @_;

  # Initialize status tracking
  my $exeFail  = 0;
  my $died     = 0;
  my $core     = 0;
  my $exitCode = 0;

  # Run command
  if(!defined($returnVal)) {
    system($command);
  } else {
    $$returnVal = `$command`;
  }

  # Check status
  if ($? == -1) {
    $exeFail = 1;
  } elsif ($? & 127) {
   $died = ($? & 127);
   $core = 1 if($? & 128);
  } else {
    $exitCode = $? >> 8;
  }

  my $problem = 0;
  if($exeFail || $died || $exitCode) {
    print STDERR "$0: problem encountered running command \"$command\"\n";
    if($exeFail) {
      print STDERR "Failed to execute command: $!\n";
    } elsif ($died) {
      print STDERR sprintf("Child died with signal %d, %s coredump\n", $died,  $core ? 'with' : 'without');
    } else {
      print STDERR "Child exited with value $exitCode\n";
    }
    $problem = 1;
  }

  return($problem);
}













