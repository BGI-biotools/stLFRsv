use strict;
use warnings;

die "usage: $0 <lnd.all file list> <bin_size> <extend_bin_num> <out_file>\n" unless @ARGV == 4;

my ($file,$bin,$ext,$out)=@ARGV;
my @name;
open IN,"$file" or die $!;
while(<IN>){
	chomp;
	push @name,$_;
}
close IN;

my %con;
foreach my $c(@name){
	open IN,"$c" or die $!;
	while(<IN>){
		chomp;
		my @t=split;
		my ($chr1,$pos1,$ori1,$chr2,$pos2,$ori2)=@t[3..8];
		if($chr1 eq $chr2){
			if($pos1<$pos2){
				my $ne="$ori1$ori2";
				unless (exists $con{"$chr1-$chr2"}{$pos1}{$pos2}){
					$con{"$chr1-$chr2"}{$pos1}{$pos2}{"RL"}=0;
					$con{"$chr1-$chr2"}{$pos1}{$pos2}{"LR"}=0;
					$con{"$chr1-$chr2"}{$pos1}{$pos2}{"LL"}=0;
					$con{"$chr1-$chr2"}{$pos1}{$pos2}{"RR"}=0;
				}
				$con{"$chr1-$chr2"}{$pos1}{$pos2}{$ne}++;
			}else{
				my $ne="$ori2$ori1";
				unless (exists $con{"$chr2-$chr1"}{$pos2}{$pos1}){
					$con{"$chr2-$chr1"}{$pos2}{$pos1}{"RL"}=0;
					$con{"$chr2-$chr1"}{$pos2}{$pos1}{"LR"}=0;
					$con{"$chr2-$chr1"}{$pos2}{$pos1}{"LL"}=0;
					$con{"$chr2-$chr1"}{$pos2}{$pos1}{"RR"}=0;
				}
				$con{"$chr2-$chr1"}{$pos2}{$pos1}{$ne}++;
			}
		
		}else{
			my $chra=$chr1;
			my $chrb=$chr2;
			$chra=~ s/^chr//;
			$chrb=~ s/^chr//;
			$chra=23 if $chra eq "X";
			$chra=24 if $chra eq "Y";
			$chra=25 if $chra eq "M";
			$chrb=23 if $chrb eq "X";
			$chrb=24 if $chrb eq "Y";
			$chrb=25 if $chrb eq "M";
			if($chra=~ /^\d+$/ and $chrb=~ /^\d+$/){
				if($chra < $chrb){
					my $ne="$ori1$ori2";
					unless (exists $con{"$chr1-$chr2"}{$pos1}{$pos2}){
						$con{"$chr1-$chr2"}{$pos1}{$pos2}{"RL"}=0;
						$con{"$chr1-$chr2"}{$pos1}{$pos2}{"LR"}=0;
						$con{"$chr1-$chr2"}{$pos1}{$pos2}{"LL"}=0;
						$con{"$chr1-$chr2"}{$pos1}{$pos2}{"RR"}=0;
					}
					$con{"$chr1-$chr2"}{$pos1}{$pos2}{$ne}++;
				
				}else{
					my $ne="$ori2$ori1";
					unless (exists $con{"$chr2-$chr1"}{$pos2}{$pos1}){
						$con{"$chr2-$chr1"}{$pos2}{$pos1}{"RL"}=0;
						$con{"$chr2-$chr1"}{$pos2}{$pos1}{"LR"}=0;
						$con{"$chr2-$chr1"}{$pos2}{$pos1}{"LL"}=0;
						$con{"$chr2-$chr1"}{$pos2}{$pos1}{"RR"}=0;
					}
					$con{"$chr2-$chr1"}{$pos2}{$pos1}{$ne}++;
				
				}
			}else{
				if($chra lt $chrb){
					my $ne="$ori1$ori2";
					unless (exists $con{"$chr1-$chr2"}{$pos1}{$pos2}){
						$con{"$chr1-$chr2"}{$pos1}{$pos2}{"RL"}=0;
						$con{"$chr1-$chr2"}{$pos1}{$pos2}{"LR"}=0;
						$con{"$chr1-$chr2"}{$pos1}{$pos2}{"LL"}=0;
						$con{"$chr1-$chr2"}{$pos1}{$pos2}{"RR"}=0;
					}
					$con{"$chr1-$chr2"}{$pos1}{$pos2}{$ne}++;
				}else{
					my $ne="$ori2$ori1";
					unless (exists $con{"$chr2-$chr1"}{$pos2}{$pos1}){
						$con{"$chr2-$chr1"}{$pos2}{$pos1}{"RL"}=0;
						$con{"$chr2-$chr1"}{$pos2}{$pos1}{"LR"}=0;
						$con{"$chr2-$chr1"}{$pos2}{$pos1}{"LL"}=0;
						$con{"$chr2-$chr1"}{$pos2}{$pos1}{"RR"}=0;
					}
					$con{"$chr2-$chr1"}{$pos2}{$pos1}{$ne}++;
				}
			}	
		}
	}
	close IN;
}

open OUT,">$out" or die $!;
my @result;
foreach my $key(keys %con){
	my ($chr1,$chr2)=split(/-/,$key);
	
	my @t_pos=sort {$a <=> $b} keys %{$con{$key}};
	while(@t_pos >0){
		my $c_pos=shift @t_pos;
		my @t_pos2=sort {$a <=> $b} keys %{$con{$key}{$c_pos}};
		if(@t_pos2 > 0){
			my @bat;
			my $count=0;
			my $m_pos=shift @t_pos2;
			push @bat,[($c_pos,$m_pos,$con{$key}{$c_pos}{$m_pos}{"RL"},$con{$key}{$c_pos}{$m_pos}{"LR"},$con{$key}{$c_pos}{$m_pos}{"LL"},$con{$key}{$c_pos}{$m_pos}{"RR"})];
			delete $con{$key}{$c_pos}{$m_pos};
			while($count != @bat){
				$count=@bat;
				&search(\@bat,$key);	
			}
			
			my (@pos1,@pos2,@num);
			@num=(0,0,0,0);
			foreach my $r (@bat){
				push @pos1,$r->[0];
				push @pos2,$r->[1];
				$num[0]+=$r->[2];
				$num[1]+=$r->[3];
				$num[2]+=$r->[4];
				$num[3]+=$r->[5];
			}
			
			@pos1=sort {$a <=> $b} @pos1;
			@pos2=sort {$a <=> $b} @pos2;
			print OUT "$chr1\t$pos1[0]\t$pos1[-1]\t$chr2\t$pos2[0]\t$pos2[-1]\t",join("\t",@num),"\n";	
		}else{
			delete $con{$key}{$c_pos};
		}
		@t_pos=sort {$a <=> $b} keys %{$con{$key}};
	}
	
	# foreach my $p1(sort {$a <=> $b} keys %{$con{$key}}){
		# foreach my $p2(sort {$a <=> $b} keys%{$con{$key}{$p1}}){
			# my @num;
			# $num[0]=$con{$key}{$p1}{$p2}{"RL"};
			# $num[1]=$con{$key}{$p1}{$p2}{"LR"};
			# $num[2]=$con{$key}{$p1}{$p2}{"LL"};
			# $num[3]=$con{$key}{$p1}{$p2}{"RR"};
			
			# print "$chr1\t$p1\t$chr2\t$p2\t".join("\t",@num)."\n";
		# }
	# }
}

close OUT;

sub search{
	my $ref=$_[0];
	my $key=$_[1];
	my (@pos1,@pos2);
	
	foreach my $r(@{$ref}){
		push @pos1,$r->[0];
		push @pos2,$r->[1];
	}
	
	@pos1=sort {$a <=> $b} @pos1;
	@pos2=sort {$a <=> $b} @pos2;
	
	my $s1=$pos1[0] - $bin*$ext;
	my $e1=$pos1[-1] + $bin*$ext;
	my $s2=$pos2[0] - $bin*$ext;
	my $e2=$pos2[-1] + $bin*$ext;
	
	for(my $i=$s1;$i<=$e1;$i+=$bin){
		for (my $j=$s2;$j<=$e2;$j+=$bin){
			if(exists $con{$key}{$i}{$j}){
				push @{$ref},[($i,$j,$con{$key}{$i}{$j}{"RL"},$con{$key}{$i}{$j}{"LR"},$con{$key}{$i}{$j}{"LL"},$con{$key}{$i}{$j}{"RR"})];
				delete $con{$key}{$i}{$j};
			}
		}
	}
}









