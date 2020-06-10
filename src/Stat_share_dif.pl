#! /usr /bin/perl -w
use strict;
use Cwd;
use FindBin qw($Bin);
die "Usage: perl $0 bam mapq bin chr1 start1 end1 chr2 start2 end2 flank prefix withchr cut type\n" if @ARGV !=14;

my ($bam,$mapq,$bin,$chr1,$start1,$end1,$chr2,$start2,$end2,$flank,$pre,$withchr,$cut,$type)=@ARGV;
my (%hash1,%hash2);
my $path = $Bin;
$cut=1-$cut;
if($withchr!=0){
	$chr1="chr".$chr1;
	$chr2="chr".$chr2;
}
my($s1,$e1,$s2,$e2);
$s1=$start1-$flank;
if($s1<0){
	$s1=0;
}
$e1=$end1+$flank;
$s2=$start2-$flank;
if($s2<0){
	$s2=0;
}
$e2=$end2+$flank;
open IN,"samtools view -F 1024 -q $mapq $bam $chr1:$s1-$e1|"or die $!;
while(<IN>){
	chomp;
	my @array=split;
	my $bar;
    if($type==0){
        $bar=(split/#/,$array[0])[-1];
        next if ($bar eq "0_0_0");
    }else{
         $_=~/BX:Z:(.*)-1/;
         $bar=$1;
    }

	my $num=int($array[3]/$bin)*$bin;
	if(exists $hash1{$num}){
		push @{$hash1{$num}},$bar;
	}else{
		$hash1{$num}=[];
		push @{$hash1{$num}},$bar;
	}
}
close IN;

open IN2,"samtools view -F 1024 -q $mapq $bam $chr2:$s2-$e2|"or die $!;
while(<IN2>){
	chomp;
	my @array=split;
	my $bar;
    if($type==0){
        $bar=(split/#/,$array[0])[-1];
        next if ($bar eq "0_0_0");
    }else{
         $_=~/BX:Z:(.*)-1/;
         $bar=$1;
    }
	my $num=int($array[3]/$bin)*$bin;
	if(exists $hash2{$num}){
		push @{$hash2{$num}},$bar;
	}else{
		$hash2{$num}=[];
		push @{$hash2{$num}},$bar;
	}
}
close IN2;

my %share;

for my $k ( sort {$a <=>$b} keys %hash1){
	my %uniq;
	@uniq{@{$hash1{$k}}} = ( );
	my @uniq_array=sort keys %uniq;
	$hash1{$k}=\@uniq_array;
}


for my $j ( sort {$a <=>$b} keys %hash2){
	my %uniq;
	@uniq{@{$hash2{$j}}} = ( );
	my @uniq_array=sort keys %uniq;
	$hash2{$j}=\@uniq_array;
}

for my $kk ( sort {$a <=>$b} keys %hash1){
	for my $jj ( sort {$a <=>$b} keys %hash2){
		for my $bar( @{$hash2{$jj}}){
			if ((grep m/$bar/,@{$hash1{$kk}})>0){
				$share{"$kk.$jj"}+=1;
			}
		}
	}
}


open OUT,">$pre.stat" or die $!;
for my $n( keys %share){
	my @info=split/\./,$n;
	print OUT "$info[0]\t$info[1]\t$share{$n}\n";
}
close OUT;
my $title=$chr1.":".$start1."-".$end1."--".$chr2.":".$start2."-".$end2;
`Rscript $path/plot_heatmap_dif.R $pre.stat $pre $title $s1 $e1 $s2 $e2 $start1 $start2 $cut $bin`;

