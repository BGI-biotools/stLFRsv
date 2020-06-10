use warnings;
use strict;
use Bio::DB::HTS;

die "$0 <refrence file> <bam file> <phase info> <out barcode info>\n" if @ARGV !=4;

my $ref=$ARGV[0];
my $bam=$ARGV[1];
my $phase=$ARGV[2];

my $sam = Bio::DB::HTS->new(-bam  => $bam,-fasta=> $ref,);
my %hash;
open IN,"$phase" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	my $seq_id=$t[0];
	my $pos=$t[1];
	my $ref=$t[2];
	my $alt=$t[3];
	my ($PAT,$MAT)=split(/\|/,$t[4]);
	my $cindel=length($alt)-length($ref);
	my $callback = sub {
        my ($seqid,$ppos,$pileup) = @_;
		return if $ppos !=$pos;
		foreach my $p(@$pileup){
			my $indel=$p->indel;
			my $a=$p->alignment;
			my $flag=$a->flag;
			next if $flag & 0x400;
			#my $barcode=$a->aux_get("BX");
			my $name=$a->qname;
			my $barcode=(split(/#/,$name))[-1];
			next if $barcode eq "0_0_0";
			if(! exists $hash{$barcode} ){
				$hash{$barcode}{"PAT"}=0;
				$hash{$barcode}{"MAT"}=0;
			}
			my $re;
			if($cindel == 0){##snp
				next unless $indel == 0;
				my $base=substr($a->qseq,$p->qpos,1);
				$re= $base eq $ref ? 0:1;
			}else{ ##indel
				if($indel!= 0){
					next unless $indel == $cindel;
					$re=1;
				}else{
					$re=0;
				}
			}
			if($re == $PAT){
				$hash{$barcode}{"PAT"}++;
			}else{
				$hash{$barcode}{"MAT"}++;
			}
		}
	};
	$sam->fast_pileup("$seq_id:$pos-$pos",$callback);
}
close IN;

open OUT,">$ARGV[3]" or die $!;
foreach my $bar(sort {$a cmp $b} keys %hash){
	my $sum=$hash{$bar}{"PAT"}+$hash{$bar}{"MAT"};
	next if $sum == 0;
	my $ratio=$hash{$bar}{"PAT"}/$sum;
	my $mark;
	if($ratio ==1){
		$mark="PAT";
	}elsif($ratio == 0){
		$mark ="MAT";
	}else{
		if($sum >=3){
			if($ratio >= 0.9){
				$mark="PAT";
			}elsif($ratio <=0.1){
				$mark="MAT";
			}else{
				$mark="NAN";
			}
		}else{
			$mark="NAN";
		}
	}
	next if $mark eq "NAN";
	print OUT "$bar\t",$hash{$bar}{"PAT"},"\t",$hash{$bar}{"MAT"},"\t",$mark,"\n";
}
close OUT;
