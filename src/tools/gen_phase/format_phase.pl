use warnings;
use strict;

die "Usage: $0 hapcut_out_phase out_merged_region out_var_phase\n" unless @ARGV == 3;
$|=1;
my $infile=$ARGV[0];
my $out_region=$ARGV[1];
my $out_var=$ARGV[2];

open IN,"$infile" or die $!;
open OR,">$out_region" or die $!;
open OV,">$out_var" or die $!;

# my @last;
my @cur;
my @check;
while(<IN>){
	chomp;
	if($_=~ /^\*/){
		if(@cur > 1){
			if(@check >0){
				my $end=0;
				foreach my $ref(@check){
					my $tmp=(split(/\t/,$ref->[-1]))[4];
					$end=$tmp if $tmp > $end;
				}
				
				my $start=(split(/\t/,$cur[1]))[4];
				if($end <= $start){
					&out_region(\@check);
					@check=();
				}
			}
			push @check,[@cur];
			# @last=@cur;
		}
		@cur=();
	}else{
		push @cur,$_ unless $_=~ /-/;
	}
}

if(@cur > 1){
	if(@check >0){
		my $end=0;
		foreach my $ref(@check){
			my $tmp=(split(/\t/,$ref->[-1]))[4];
			$end=$tmp if $tmp > $end;
		}
		my $start=(split(/\t/,$cur[1]))[4];
		if($end <= $start){
			&out_region(\@check);
			@check=();
		}
	}
	push @check,[@cur];
}
if(@check >0){
	&out_region(\@check);
}
close IN;
close OR;
close OV;

sub out_region{
	my @data=@{$_[0]};
	if(@data == 1){
		my $start=(split(/\t/,$data[0]->[1]))[4];
		my $end=(split(/\t/,$data[0]->[-1]))[4];
		# @ all output
		print OR "$start\t$end\n";
		shift @{$data[0]};
		foreach my $line(@{$data[0]}){
			my @tmp=split(/\t/,$line);
			print OV "$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[1]|$tmp[2]\n";
		}
	}else{
		my %index;
		my @count;
		my $total;
		for(my $i=0;$i<=$#data;$i++){
			$count[$i]=@{$data[$i]};
			$total+=$#{$data[$i]};	
		}
		
		for(my $i=0;$i<=$#data;$i++){
			my $ratio=$count[$i]/$total;
			next if $count[$i] < 3;
			next if $ratio < 0.1;
			shift @{$data[$i]};# the BLOCK line
			my $spos=(split(/\t/,$data[$i]->[0]))[4];
			$index{$i}=$spos;
		}
		
		while((keys %index) > 0){
			if((keys %index) == 1){
				my $in=(keys %index)[0];
				my $start=(split(/\t/,$data[$in]->[0]))[4];
				my $end=(split(/\t/,$data[$in]->[-1]))[4];
				# @ all output
				print OR "$start\t$end\n";
				foreach my $line(@{$data[$in]}){
					my @tmp=split(/\t/,$line);
					print OV "$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[1]|$tmp[2]\n";
				}
				delete $index{$in};
			}else{
				my @tmppos=sort{
					$index{$a} <=> $index{$b}
				}(keys %index);
				my $min=$tmppos[0];
				my $max=$tmppos[1];
				
				my @record;
				my $spos=$index{$min};
				while(($spos < $index{$max}) and (@{$data[$min]} > 0)){
					push @record, shift @{$data[$min]};
					$spos=(split(/\t/,$data[$min]->[0]))[4] if (@{$data[$min]} > 0);
				}
				if(@record >= 2){
					my $start=(split(/\t/,$record[0]))[4];
					my $end=(split(/\t/,$record[-1]))[4];
					print OR "$start\t$end\n";
					foreach my $line(@record){
						my @tmp=split(/\t/,$line);
						print OV "$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[1]|$tmp[2]\n";
					}
				}
				#remamin at least 3 element
				if(@{$data[$min]} > 1){
					$index{$min}=$spos;
				}else{
					delete $index{$min};
				}	
			}
		}
	}
}