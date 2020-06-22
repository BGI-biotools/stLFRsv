#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);
die "perl $0 stLFRsv.final smoove.vcf flank ratio break sample filt" if @ARGV<7;
my $sv=shift;
my $smoove=shift;
my $flank=shift;
my $ratio=shift;
my $break=shift;
my $sample=shift;
my $filt=shift;
$filt||=1;
my $path=$Bin;
open IN,"<$sv" or die $!;

<IN>;
my %hash;
my %type;
while(<IN>){
	chomp;
	my @l=split;
    if($filt==1){
		next unless $l[11] eq "PASS";
    }elsif($filt==2){
		next unless ($l[11] eq "PASS" || $l[11] eq "PASS|COMMON");
        if($l[11] eq "PASS|COMMON"){
        	$l[11]="PASS";
		}
	}
	my $gt=(split/:/,$l[13],)[1];
    if($gt=~/[A-Z]/){
    	$gt="./.";
	}
	$l[4]=~s/chr//;
	$l[6]=~s/chr//;
	if($l[10] eq "DEL"|| $l[10] eq "DUP"||$l[10]=~/INV/){
		$l[10]=~/(\w+)/;
		my $type=$1;
		my $len=$l[7]-$l[5];
		next if ($len<=$break);
		if($type eq "DEL"){
			$len="-".$len;
		}elsif($type eq "DUP"){
			$len="+".$len;
		}
		push @{$hash{$l[4]}{$l[5]}},"$l[4]\t$l[5]\t.\tN\t<$type>\t.\t$l[11]\tSVTYPE=$type;SVLEN=$len;END=$l[7];SOURCE=stLFRsv\tGT\t$gt\n";
		push @{$type{"$type\_$l[4]"}},[$l[5],$l[7],abs($len),$l[11]];
	}else{
		my ($type1,$type2);
		if($l[10] eq "TRA1"){
			$type1="N\[$l[6]:$l[7]\[";
			$type2="\]$l[4]:$l[5]\]N";
		}elsif($l[10] eq "TRA2"){
			$type2="N\[$l[4]:$l[5]\[";
			$type1="\]$l[6]:$l[7]\]N";
		}elsif($l[10] eq "TRA3"){
			$type1="\[$l[6]:$l[7]\[N";
			$type2="N\]$l[4]:$l[5]\]";
		}elsif($l[10] eq "TRA4"){
			$type2="\[$l[4]:$l[5]\[N";
			$type1="N\]$l[6]:$l[7]\]";
		}
		push @{$hash{$l[4]}{$l[5]}},"$l[4]\t$l[5]\t.\tN\t$type1\t.\t$l[11]\tSVTYPE=BND;SOURCE=stLFRsv\tGT\t$gt\n";	
		push @{$hash{$l[6]}{$l[7]}},"$l[6]\t$l[7]\t.\tN\t$type2\t.\t$l[11]\tSVTYPE=BND;SOURCE=stLFRsv\tGT\t$gt\n";	
		my($c1,$c2);
		$c1=$l[9];
		if($c1 eq "RL"){
			$c2="LR";
		}elsif($c1 eq "LR"){
			$c2="RL";
		}else{
			$c2=$c1;
		}
		push @{$type{"BND_$l[4]"}},[$l[5],$l[6],$l[7],$c1,$l[11]];
		push @{$type{"BND_$l[6]"}},[$l[7],$l[4],$l[5],$c2,$l[11]];
	}

}
close IN;

if($smoove=~/\.gz$/){
	open IN2,"gzip -dc $smoove| " or die $!;
}else{
	open IN2,"<$smoove" or die $!;
}

while(<IN2>){
	next if /^#/;
	chomp;
	my @l=split;
	$_=~/SVTYPE=(.*?);/;
	my $type=$1;
	my $merge=0;
	my $info;
	$l[0]=~s/chr//;
	if($type=~/DEL|DUP|INV/){
		$_=~/SVLEN=(.*?);END=(.*?);/;
		my $len=abs($1);
		next if $len>$break;
		next unless $len>50;
		my $e=$2;
		my $s=$l[1];
		if(exists $type{"$type\_$l[0]"}){
			my($tmp_s,$tmp_e);
			for my $k(0..scalar(@{$type{"$type\_$l[0]"}})-1){
				$tmp_s=$type{"$type\_$l[0]"}[$k][0];
				$tmp_e=$type{"$type\_$l[0]"}[$k][1];
				my $maxs=&max($s,$tmp_s);
				my $mine=&min($e,$tmp_e);
				my $overlap=$mine-$maxs;
				next if $overlap<=0;
				if($overlap/$len >=$ratio && $overlap/$type{"$type\_$l[0]"}[$k][2] >=$ratio){
					$l[6]=$type{"$type\_$l[0]"}[$k][-1];
					$l[7]=~s/STRANDS=.*//;
					$l[7].="SOURCE=merge";
					$l[8]=~s/:GQ.*//;
					$l[9]=(split(":",$l[9]))[0];
					my $info=join"\t",@l;
					push @{$hash{$l[0]}{$s}},"$info\n";
					delete $hash{$l[0]}{$tmp_s};	
					$merge=1;
					last;
				}
			}
		}
		if($merge==0){
			$l[7]=~s/STRANDS=.*//;
			$l[7].="SOURCE=smoove";
			$l[8]=~s/:GQ.*//;
			$l[9]=(split(":",$l[9]))[0];
			$info=join"\t",@l;
			push @{$hash{$l[0]}{$s}},"$info\n";
		}
	}elsif($type eq "BND"){
		$l[4]=~/(N\[|N\]|\[N|\]N)/;
		my $c;
		if($1 eq "N["){
			$c="RL";
		}elsif($1 eq "N]"){
			$c="RR";
		}elsif($1 eq "[N"){
			$c="LL";
		}elsif($1 eq "]N"){
			$c="LR";
		}
		$l[4]=~s/chr//;
		if(exists $type{"$type\_$l[0]"}){
			for my $k(0..scalar(@{$type{"$type\_$l[0]"}})-1){
				if($type{"$type\_$l[0]"}[$k][3] eq $c && abs($l[1] - $type{"$type\_$l[0]"}[$k][0])<=$flank ){
					$l[4]=~/(\]|\[)(.*):(.*)(\]|\[)/;
					my $chr=$2;
					my $pos=$3;
					if($chr eq $type{"$type\_$l[0]"}[$k][1] && abs($pos - $type{"$type\_$l[0]"}[$k][2])<=$flank){
						$l[6]=$type{"$type\_$l[0]"}[$k][-1];
						$l[7]=~s/STRANDS=.*//;
						$l[7].="SOURCE=merge";
						$l[8]=~s/:GQ.*//;
						$l[9]=(split(":",$l[9]))[0];
						my $info=join"\t",@l;
						push @{$hash{$l[0]}{$l[1]}},"$info\n";						
						$merge=1;
						last;
					}
				}
			 	
			}
		}	
		if($merge==0){
			$l[7]=~s/STRANDS=.*//;
			$l[7].="SOURCE=smoove";
			$l[8]=~s/:GQ.*//;
			$l[9]=(split(":",$l[9]))[0];
			$info=join"\t",@l;
			push @{$hash{$l[0]}{$l[1]}},"$info\n";
		}
	}else{
		$l[7]=~s/STRANDS=.*//;
		$l[7].="SOURCE=smoove";
		$l[8]=~s/:GQ.*//;
		$l[9]=(split(":",$l[9]))[0];
		$info=join"\t",@l;
		push @{$hash{$l[0]}{$l[1]}},"$info\n";
	}
}
close IN2;

open OUT,">merge.dels.vcf" or die $!;
open OUT2,">merge.other.sv.vcf" or die $!;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $y=$year+1900;
my $m=$mon+1;
if(length($m)==1){
    $m="0".$m;
}
my $d=$mday;
if(length($d)==1){
    $d="0".$d;
}
my $date=$y.$m.$d;
open HE,"<$path/merge.dels.header" or die $!;

while(<HE>){
    if($_=~/fileDate/){
        $_=~s/=(\d+)/=$date/;
    }elsif($_!~/^##/){
        chomp;
        $_.="\t$sample\n";
    }
    print OUT $_;
}
close HE;

open HE2,"<$path/merge.others.header" or die $!;

while(<HE2>){
    if($_=~/fileDate/){
        $_=~s/=(\d+)/=$date/;
    }elsif($_!~/^##/){
        chomp;
        $_.="\t$sample\n";
    }
    print OUT2 $_;
}
close HE2;

my ($c1,$c2);
my $line;
$c1=$c2=0;
for my $k(1..22){
	for my $kk (sort {$a<=>$b }keys %{$hash{$k}}){
		for my $kkk (@{$hash{$k}{$kk}}){
			if($kkk=~/DEL/){
				$c1+=1;
				chomp $kkk;
				my @t=split/\t/,$kkk;
				$t[2]=$c1;
				$line=join"\t",@t;
				print OUT "$line\n";
			}else{
				$c2+=1;
				chomp $kkk;
				my @t=split/\t/,$kkk;
				$t[2]=$c2;
				$line=join"\t",@t;
				print OUT2 "$line\n";
			}
		}
	}
}
for my $k (sort {$a<=>$b} keys %{$hash{"X"}}){
	for my $kk (@{$hash{"X"}{$k}}){
		if($kk=~/DEL/){
			$c1+=1;
			chomp $kk;
			my @t=split/\t/,$kk;
			$t[2]=$c1;
			$line=join"\t",@t;
			print OUT "$line\n";
		}else{
			$c2+=1;
			chomp $kk;
			my @t=split/\t/,$kk;
			$t[2]=$c2;
			$line=join"\t",@t;
			print OUT2 "$line\n";
		}
	}
}

for my $k (sort {$a<=>$b} keys %{$hash{"Y"}}){
	for my $kk (@{$hash{"Y"}{$k}}){
		if($kk=~/DEL/){
			$c1+=1;
			chomp $kk;
			my @t=split/\t/,$kk;
			$t[2]=$c1;
			$line=join"\t",@t;
			print OUT "$line\n";
		}else{
			$c2+=1;
			chomp $kk;
			my @t=split/\t/,$kk;
			$t[2]=$c2;
			$line=join"\t",@t;
			print OUT2 "$line\n";
		}
	}

}

sub min {
	my $min = shift;
	for (@_) {
		$min = $_ if $_ < $min;
	}
	return $min;
}


sub max {
	my $max = shift;
	for (@_) {
		$max = $_ if $_ > $max;
	}
	return $max;
}
