#! /usr /bin/perl -w
use strict;
use Cwd;
use FindBin qw($Bin);

die "Usage: perl $0 bam mapq bin chr start end flank prefix withchr cut type \n" if @ARGV != 11;

my ($bam,$mapq,$bin,$chr,$start,$end,$flank,$pre,$withchr,$cut,$type)=@ARGV;
my %hash;
if($withchr!=0){
	$chr="chr".$chr;
}
$cut=1-$cut;
my $path = $Bin;
my($s,$e);
$s=$start-$flank;
if($s<0){
	$s=0;
}
$e=$end+$flank;
open IN,"samtools view -F 1024 -q $mapq $bam $chr:$s-$e|"or die $!;
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
	if(exists $hash{$num}){
		push @{$hash{$num}},$bar;
	}else{
		$hash{$num}=[];
		push @{$hash{$num}},$bar;
	}
}
close IN;
my %share;

for my $k ( sort {$a <=>$b} keys %hash){
	my %uniq;
	@uniq{@{$hash{$k}}} = ( );
	my @uniq_array=sort keys %uniq;
	$hash{$k}=\@uniq_array;
}

for my $k ( sort {$a <=>$b} keys %hash){
	for my $j ( sort {$a <=>$b} keys %hash){
		for my $bar( @{$hash{$j}}){
			if ((grep m/$bar/,@{$hash{$k}})>0){
				$share{"$k.$j"}+=1;
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
my $title=$chr.":".$start."-".$end;
`Rscript $path/plot_heatmap.R $pre.stat $pre $title $s $e $start $end $cut $bin`;


