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
my $sd;
my $p_th;
my $phase_dir;
my $black_list;
my $control_list;
my $id_num;
my $Smerge;
my $ncpu;
my $human;
my $mergemax;
my $qc1;
my $qc2;
my $qc3;
my $rlen;
my $mlen;
my $sp;
my $cn;
my $bin_r;
my $help;


unshift @PATH, "$path/bin";
my $Basename=basename($0);
my $USAGE = qq{
Name:
	$Basename
	version 2.2
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
	-gap <int> define the gap size which should be considered as different segment.
	-size <int> output SV length.[default 20000]
	-is <int> proper IS size for read pair library, read pairs with too large IS will be abandoned.[default 300]
	-bin <int> bin size for cluster the segments.
	-merge1 <int> N continue outline bins could be considered as the same break point and will be merged into one evidence.
	-merge2 <int> SVs nearby under N binsize will be considered as one event.[default 5]
	-mmax <int> the max SVs allowed in one event.[default 4]
	-low <int> lowest shared barcode counts threshold.[default 4]
	-sd <int> break ends with a depth higher than avg_dep+N*sd will be considered as candidates.[default 3]
	-p_th <float> break ends significantly high with P value lower than this threshold will be considered as candidates.[default 0.1]
	-phase <string> formatted phase result directory including phased barcode and region by chromosome.[default NULL]
	-bl <string> black list file(BED format).[default NULL]
	-cl <string> sorted control list file(BEDPE format).[default NULL](Be sure the chromosome and position are sorted in one line!!!)
	-sc <int> allow max sv counts for the same position in one direction.[default 4]
	-human <Y/N> for Homo sapiens,keep only [1234567890XYM] chromosome.[default N]
	-qc1 <float> valid read pair ratio for SV detection.[default 0.60]
	-qc2 <float> average read pair count for one barcode.[default 30]
	-qc3 <float> average segment end count for one bin.[default 15]
	-sp <float> sample percentage for DNA fragment length statistic.[default 0.2]
	-cn <int> sample count for read pair distance statistic.[default 20000000]
	-rlen <int> read length of one read.[default 100]
	-mlen <int> physical limit for the long DNA segment.[default 400000]
	-help Show this message.
};

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
"low=i" =>\$low,
"sd=i" =>\$sd,
"p_th=f" =>\$p_th,
"mmax=i" =>\$mergemax,
"phase=s" => \$phase_dir,
"bl=s" => \$black_list,
"cl=s" => \$control_list,
"sc=i" => \$id_num,
"human=s" => \$human,
"qc1=f" => \$qc1,
"qc2=f" => \$qc2,
"qc3=f" => \$qc3,
"sp=f" => \$sp,
"cn=i" => \$cn,
"rlen=i" => \$rlen,
"mlen=i" => \$mlen
);

$ncpu ||= 1;
$tmp ||= "/tmp";
$bar_th||=8;
$seg_th ||= 4;
$size ||= 20000;
$is ||=300;
$low ||=4;
$sd ||=3;
$p_th ||=0.1;
$phase_dir ||="NULL";
$black_list||="NULL";
$control_list||="NULL";
$id_num||=4;
$Smerge||=5;
$human||="N";
$mergemax ||=4;
$qc1 ||=0.6;
$qc2 ||=30;
$qc3 ||=15;
$rlen ||=100;
$mlen ||=400000;
$sp ||=0.2;
$cn ||=20000000;

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

my $stat=0;
if(defined $gap and defined $bin and defined $Nmerge){
	print STDERR "Use the user specific parameters\n";
	print STDERR "Bin size: $bin\n";
	print STDERR "Merge windows: $Nmerge\n";
	print STDERR "Gap size: $gap\n";
}else{
	$stat=1;
}

if($stat){
	print STDERR "User didn't specify the bin, merge1 and gap\n";
	print STDERR "Will set these parameters automatically\n";
	print STDERR "Gather read pair distances form the largest contig\n";
	$line="barcode-sort $bam $ncpu $out $bamname $bar_th $seg_th $mlen $is $human Y $cn";
	if(executeSystemCall($line)){
		die "Failure during creating the sorted barcoded file $out/$bamname.gap\n";
	}
	
	my $R = Statistics::R->new();
	my $cmd=qq{
a = read.table("$out/$bamname.gap")
x<-a[,1]
x<-x[which(x<$mlen)]
x<-x[which(x>$rlen)]
len1<-quantile(x,0.65)
len2<-quantile(x,0.93)
len3<-quantile(x,0.98)
};
	my ($s1,$s2,$s3);
	$R->run($cmd);
	$s1=$R->get('len1');
	$s2=$R->get('len2');
	$s3=$R->get('len3');
	
	$Nmerge=int($s2/$s1+0.5);
	$bin=int($s1/100+0.5)*100;
	$gap=int($s3/100+0.5)*100;
	print STDERR "According to the read pair distance statistic:\n";
	print STDERR "bin size is set to $bin\n";
	print STDERR "merge bin number is set to $Nmerge\n";
	print STDERR "gap size is set to $gap\n";
	
	$bin_r=0.65;
}

print STDERR "Creat the sorted barcoded file\n";
$line="barcode-sort $bam $ncpu $out $bamname $bar_th $seg_th $gap $is $human N $cn";
if(executeSystemCall($line)){
	die "Failure during creating the sorted barcoded file $out/$bamname.sbf\n";
}

open IN,"$out/$bamname.stat" or die $!;
my $tl=<IN>;
close IN;
chomp $tl;
my ($tread,$vread,$vbar,$vseg)=split(/\t/,$tl);
my $q1=$vread/$tread;
my $q2=$vread/$vbar;
my $q3=$vseg*2/($reflen/$bin);
if($q1 < $qc1){
	print STDERR "Warning: valid read pair ratio $q1 is lower than $qc1, the result may be unreliable!\n";
}

if($q2 < $qc2){
	print STDERR "Warning: average read pair on one barcode $q2 is lower than $qc2, the result may be unreliable!\n";
}

if($q3 < $qc3){
	print STDERR "Warning: average segment end in one bin $q3 is lower than $qc3, the result may be unreliable!\n";
}


unless($stat){
	my $R = Statistics::R->new();
	my $cmd=qq{
a = read.table("$out/$bamname.all.gap")
x<-a[,1]
x<-x[which(x<$mlen)]
x<-x[which(x>$rlen)]
ratio<-length(which(x < $bin))/length(x)
};
	$R->run($cmd);
	$bin_r=$R->get('ratio');
}


###########################
######freq for all#########
my $bar_max=$vbar*$sp;
my %bar_freq;
$bar_max= 2000000 if $bar_max < 2000000;
open IN,"$out/$bamname.seg" or die $!;
my @one_bar;
my $countbar=0;
while(<IN>){
	last if $countbar == $bar_max;
	chomp;
	my $clen=$_;
	if($clen == 0xFFFFFFFF){
		if(@one_bar > 0){
			&process_bar;
			@one_bar=();
			$countbar++;
		}
	}else{
		push @one_bar,$clen;
	}
}


close IN;

open OUT,">$out/$bamname.freq" or die $!;
foreach my $key(sort {$a <=> $b} keys %bar_freq){
	my $f_ratio=($bar_freq{$key}->[0]/$bar_freq{$key}->[1])*$bin_r;
	print OUT "$key\t$f_ratio\n";
}
close OUT;

######freq for high quanlity#########
$countbar=0;
@one_bar=();
%bar_freq=();
my $lbar=-1;
open OUT,">$out/$bamname.HQ.seg" or die $!;
open IN,"$out/$bamname.sbf" or die $!;
binmode(IN);
my $buf;
while(read(IN,$buf,8)){
	last if $countbar == $bar_max;
	my $bar=unpack("%64Q",$buf);

	read(IN,$buf,4);
	my $index=unpack("%32L",$buf);

	read(IN,$buf,32);
	my $name=unpack("Z32",$buf);

	read(IN,$buf,4);
	my $s=unpack("%32L",$buf);

	read(IN,$buf,4);
	my $e=unpack("%32L",$buf);

	read(IN,$buf,4);
	my $sup=unpack("%32L",$buf);
	
	if($bar != $lbar){
		if(@one_bar >0){
			&process_bar;
			@one_bar=();
			$lbar=$bar;
			$countbar++;
		}
		push @one_bar,$e-$s;
		print OUT $e-$s,"\n";
	}else{
		push @one_bar,$e-$s;
		print OUT $e-$s,"\n";
	}
}
close IN;
close OUT;

open OUT,">$out/$bamname.HQ.freq" or die $!;
foreach my $key(sort {$a <=> $b} keys %bar_freq){
	my $f_ratio=($bar_freq{$key}->[0]/$bar_freq{$key}->[1])*$bin_r;
	print OUT "$key\t$f_ratio\n";
}
close OUT;

my $seg_size;
if(1){
	my $R = Statistics::R->new();
	my $cmd=qq{
a = read.table("$out/$bamname.HQ.seg")
x<-a[,1]
len1<-quantile(x,0.65)
};
	$R->run($cmd);
	$seg_size=$R->get('len1');
	$seg_size=int($seg_size/100+0.5)*100;
	print STDERR "Most high quanlity segment length is under $seg_size.(65%)\n";
}

###########################
print STDERR "Step_1 single end cluster\n";
$line="sin-cluster $out/$bamname.sbf $out/$bamname.sin $bin";

if(executeSystemCall($line)){
	die "Failure during creating the one end file $out/$bamname.sin\n";
}

print STDERR "Step_2 link the cluster to cluster\n";
my $avg_dep=int($q3+0.5);
$line="link-id $out/$bamname.sin $out/$bamname.lnd $out/$bamname.HQ.freq $low $avg_dep $sd $p_th $gap $seg_size $size $bin $Nmerge $id_num";

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
$line="judge-link $out/$bamname.lnd.all $out/$bamname.sln $bam $out/$bamname.freq $phase_dir $out/$bamname.judge $bin $gap $Nmerge $ncpu";
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
my $ext=4*$gap;
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

sub process_bar{
	@one_bar=sort {$a <=> $b} @one_bar;
	my $total_len=0;
	foreach my $l (@one_bar){
		$total_len+=$l;
	}

	for(my $i=$bin;$i<=$mlen;$i+=$bin){
		my $cur_len=0;
		foreach my $l (@one_bar){
			if($l >= $i){
				$cur_len+=$l-$i;
			}
		}
		
		my $ratio;
		if($total_len > 0){
			$ratio=$cur_len/$total_len;
		}else{
			$ratio=0;
		}
		
		if(exists $bar_freq{$i}){
			$bar_freq{$i}->[0]+=$ratio;
			$bar_freq{$i}->[1]++;
		}else{
			$bar_freq{$i}->[0]=$ratio;
			$bar_freq{$i}->[1]=1;
		}
	}
	return;
}












