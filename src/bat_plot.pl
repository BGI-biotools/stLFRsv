use strict;
use warnings;
use FindBin qw($Bin);
use Env qw(@PATH);
use threads;
use Thread::Semaphore;
use threads::shared;
unshift @PATH, $Bin;

die "Usage: $0 bam sv_final map_quality binsize extend_len cut PASS_only out_dir ncpu bam_type\n" unless @ARGV==10;

my ($bam,$file,$mapq,$bin,$ext,$cutoff,$pass,$out,$ncpu,$type)=@ARGV;

unless(-d $out){
	`mkdir -p $out`;
}

my $semaphore=new Thread::Semaphore($ncpu);

open IN,"$file" or die $!;
while(<IN>){
	next unless /^S/;
	chomp;
	my @t=split;
	if($pass){
		next unless $t[12] eq "PASS";
	}
	
	my ($chr1,$pos1,$chr2,$pos2)=@t[4..7];
	my $ltype=$t[10];
	my $line;
	my $name="$chr1-$pos1--$chr2-$pos2-$ltype";
	if($chr1 eq $chr2){
		my $s1=$pos1-$ext;
		my $e1=$pos1+$ext;
		
		$s1=0 if $s1<0;
		
		my $s2=$pos2-$ext;
		my $e2=$pos2+$ext;
		
		$s2=0 if $s2<0;
		
		if($e1 >= $s2){
			$line="Stat_share $bam $mapq $bin $chr1 $s1 $e2 0 $out/$name 0 $cutoff $type > $out/$name.log 2>&1";
		}else{
			$line="Stat_share_dif $bam $mapq $bin $chr1 $s1 $e1 $chr2 $s2 $e2 0 $out/$name 0 $cutoff $type > $out/$name.log 2>&1";
		}
	}else{
		my $s1=$pos1-$ext;
		my $e1=$pos1+$ext;
		
		$s1=0 if $s1<0;
		
		my $s2=$pos2-$ext;
		my $e2=$pos2+$ext;
		
		$s2=0 if $s2<0;
		
		$line="Stat_share_dif $bam $mapq $bin $chr1 $s1 $e1 $chr2 $s2 $e2 0 $out/$name 0 $cutoff $type > $out/$name.log 2>&1";
		
	}
	$semaphore->down();
	my $thread=threads->new(\&one_exe,$line);
	foreach my $t (threads->list(threads::joinable)){
		my $value=$t->join();
		if($value){
			exit(1);
		}
	}
}

close IN;

&waitquit;

sub one_exe{
	my $line=$_[0];
	if(executeSystemCall($line)){
		die "Failure during run [$line]\n";
	}
	
	$semaphore->up();
	return;
}

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

sub waitquit{
	my $num=0;
	while($num<$ncpu){
		$semaphore->down();
		$num++;
		foreach my $t (threads->list(threads::joinable)){
			my $value=$t->join();
			if($value){
				exit(1);
			}
		}
	}
	$semaphore->up($ncpu);#reset
}