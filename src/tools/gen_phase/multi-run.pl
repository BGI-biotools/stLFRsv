#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Env qw(@PATH);
use threads;
use Thread::Semaphore;
use threads::shared;
unshift @PATH, $Bin;

die "Usage: $0 ref bam phase_dir phase_prefix out_dir ncpu\n" unless @ARGV==6;

my ($ref,$bam,$phase_dir,$pre,$out,$ncpu)=@ARGV;

unless(-d $out){
	`mkdir -p $out`;
}
my @n;
opendir DIR,"$phase_dir" or die $!;
foreach my $name(readdir DIR){
	next unless $name=~ /^$pre/;
	push @n,$name;
}
closedir DIR;

my $semaphore=new Thread::Semaphore($ncpu);

foreach my $name(@n){
	$semaphore->down();
	$name=~ s/$pre//;
	
	my $thread=threads->new(\&sub,$name);
	foreach my $t (threads->list(threads::joinable)){
		my $value=$t->join();
		if($value){
			exit(1);
		}
	}
}
&waitquit;

sub sub{
	my $chr=$_[0];
	my $line;
	$line="format_phase $phase_dir/$pre$chr $out/$chr.region $out/$chr.vcf";
	if(executeSystemCall($line)){
		die "Failure during run [$line]\n";
	}
	$line="get_barcode_from_phase $ref $bam $out/$chr.vcf $out/$chr.barcode.phase";
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
