use strict;
use warnings;
use Statistics::R;
use List::Util;

die "Usage: $0 <cluster file> <link id file> <low depth> <exclude depth> <low depth share> <exclude depth share> <gap size> <detect size> <bin size> <N merge bins>\n" if @ARGV != 10;

my ($infile,$outfile,$low,$ex,$lows,$exs,$gap,$size,$bin,$Nmerge)=@ARGV;

my $stat=0;
if($lows <= 0){
	die "Wrong parameter for low depth share filter!\n";
}elsif($lows < 1){
	$stat=1;
}

if(1){
	my $total_len=0;
	my %barhash;
	my @cluarray;
	my @infoarray;
	open IN,"$infile.raw" or die $!;
	while(<IN>){
		chomp;
		$total_len++;
		my ($gid,$chr,$pos,$line)=(split)[0,2,3,6];
		push @infoarray,($chr,$pos);
		foreach my $sinbar(split(/,/,$line)){
			$sinbar=~ /^(\d+)-/;
			my $bar=$1;
			push @{$barhash{$bar}},$gid;
			push @{$cluarray[$gid]},$bar;
		}
	}
	close IN;
	
	open OUT,">$infile.share" or die $!;
	open OUTI,">$infile.share.dat" or die $!;
	for(my $i=0;$i<$total_len;$i++){
		my (%tmphash,$chr1,$pos1,$chr2,$pos2,$dis);
		($chr1,$pos1)=@infoarray[2*$i,2*$i+1];
		foreach my $bar ( @{$cluarray[$i]} ){
			foreach my $j(@{$barhash{$bar}}){
				next unless $i < $j;
				($chr2,$pos2)=@infoarray[2*$j,2*$j+1];
				if($chr1 eq $chr2){
					$dis=abs($pos1-$pos2);
					next if $dis < $gap;
				}
				$tmphash{$j}{$bar}++;
			}
		}
		
		if((keys %tmphash) > 0 ){
			foreach my $j (sort {$a <=> $b} keys %tmphash){
				my $count=keys %{$tmphash{$j}};
				my $sline=join(",",(keys %{$tmphash{$j}}));
				print OUT "$i\t$j\t$count\t$sline\n";
				print OUTI "$count\n";
			}
		}	
	}
	close OUT;
	close OUTI;	
}



if($stat){
	my $R = Statistics::R->new();
	my $cmd=qq{
library(MASS)
a = read.table("$infile.share.dat")
x<-a[,1]
subx<-x[which(x>(quantile(x,$exs)))]
res=fitdistr(subx, "lognormal")
th1<-qlnorm($lows,res[[1]][1],res[[1]][2],lower.tail = TRUE)
th2<-quantile(x,$lows)
};

	$R->run($cmd);
	my $low1=$R->get('th1');
	my $low2=$R->get('th2');
	$R->stop();
	print STDERR "For both ends at $lows confidence:\n";
	print STDERR "threshold1 is set to $low1\n";
	print STDERR"threshold2 is set to $low2\n";
	print STDERR "final threshold is set to ".int(($low1+$low2)/2 + 0.5)."\n";
	if(abs($low1-$low2) >= ($low1< $low2 ? $low1:$low2)/2){
		print STDERR "Warning: too large diff between threshold1 and threshold2, the result may be unreliable!!!\n";
	}
	$lows=int(($low1+$low2)/2 + 0.5);
}

$stat=0;
if($low <= 0){
	die "Wrong parameter for low depth single filter!\n";
}elsif($low < 1){
	$stat=1;
}

if($stat){
	open IN,"$infile.raw" or die $!;
	open OUT,">$infile.single.dat" or die $!;
	while(<IN>){
		chomp;
		my @t=split;
		print OUT "$t[5]\n";
	}
	close IN;
	close OUT;
	
	my $R = Statistics::R->new();
	my $cmd=qq{
library(MASS)
a = read.table("$infile.single.dat")
x<-a[,1]
subx<-x[which(x>(quantile(x,$ex)))]
res=fitdistr(subx, "lognormal")
th1<-qlnorm($low,res[[1]][1],res[[1]][2],lower.tail = TRUE)
th2<-quantile(x,$low)
};
	$R->run($cmd);
	my $low1=$R->get('th1');
	my $low2=$R->get('th2');
	$R->stop();
	print STDERR "For single end at $low confidence:\n";
	print STDERR "threshold1 is set to $low1\n";
	print STDERR "threshold2 is set to $low2\n";
	print STDERR "final threshold is set to ".int(($low1+$low2)/2 + 0.5)."\n";
	if(abs($low1-$low2) >= ($low1< $low2 ? $low1:$low2)/2){
		print STDERR "Warning: too large diff between threshold1 and threshold2, the result may be unreliable!!!\n";
	}
	$low=int(($low1+$low2)/2 + 0.5);
}

my @sv;
if(1){
	my @infoarray;
	open IN,"$infile.raw" or die $!;
	while(<IN>){
		chomp;
		my ($chr,$pos,$ori,$count)=(split)[2,3,4,5];
		push @infoarray,($chr,$pos,$ori,$count);
	}
	close IN;
	
	my @sharearray;
	my %existhash;
	open IN,"$infile.share" or die $!;
	while(<IN>){
		chomp;
		my ($i,$j,$count,$line)=split;
		$sharearray[$i]{$j}=[($count,$line)];
	}
	close IN;
	
	foreach (my $i=0;$i<=$#sharearray;$i++){
		foreach my $j(sort {$a <=> $b} keys %{$sharearray[$i]}){
			next unless $sharearray[$i]{$j}->[0] >= $lows;
			next unless ($infoarray[4*$i+3] >= $low and $infoarray[4*$j+3] >= $low);
			next if exists $existhash{$i}{$j};
			my $chra=$infoarray[4*$i+0];
			my $chrb=$infoarray[4*$j+0];
			my $posa=$infoarray[4*$i+1];
			my $posb=$infoarray[4*$j+1];
			if($chra eq $chrb){
				next unless abs($posa-$posb) >=$size;
			}
			my @result;
			&searchBest($i,$j,\@infoarray,\@sharearray,\%existhash,\@result);
			push @sv,[@result];
		}
	}
	
}

my @new_sin;
if(1){
	my %sin;
	my @new_sv;
	my $index=0;
	foreach my $ref(@sv){
		my @result=@{$ref};
		my $jcount=0;
		#chra posa oria chrb posb orib count id bar
		my (@id1,@id2);
		foreach my $k(split(/,/,$result[7])){
			$jcount++;
			my ($i,$j)=split(/-/,$k);
			push @id1,$i;
			push @id2,$j;
		}
		
		unless(exists $sin{$result[2]}{$result[0]}{$result[1]}){
			$sin{$result[2]}{$result[0]}{$result[1]}{"sid"}=$index;
			%{$sin{$result[2]}{$result[0]}{$result[1]}{"id"}}=();
			$index++;
		}
		
		foreach my $d(@id1){
			$sin{$result[2]}{$result[0]}{$result[1]}{"id"}{$d}=1;
		}

		unless(exists $sin{$result[5]}{$result[3]}{$result[4]}){
			$sin{$result[5]}{$result[3]}{$result[4]}{"sid"}=$index;
			%{$sin{$result[5]}{$result[3]}{$result[4]}{"id"}}=();
			$index++;
		}
		
		foreach my $d(@id2){
			$sin{$result[5]}{$result[3]}{$result[4]}{"id"}{$d}=1;
		}
		
		push @new_sv,[($sin{$result[2]}{$result[0]}{$result[1]}{"sid"},$sin{$result[5]}{$result[3]}{$result[4]}{"sid"},@result,$jcount)];
	}
	
	foreach my $ori(keys %sin){
		foreach my $chr (keys %{$sin{$ori}}){
			foreach my $pos (keys %{$sin{$ori}{$chr}}){
				my @temp_id=sort {$a <=> $b} keys %{$sin{$ori}{$chr}{$pos}{"id"}};
				$new_sin[$sin{$ori}{$chr}{$pos}{"sid"}]=[($chr,$pos,$ori,join(",",@temp_id))];
			}
		}
	}	
	
	@sv=@new_sv;
}

if(1){
	my @infoarray;
	open IN,"$infile.raw" or die $!;
	while(<IN>){
		chomp;
		my @t=split;
		push @infoarray,$t[6];
	}
	close IN;
	
	open OUT,">$infile" or die $!;
	for (my $i=0;$i<=$#new_sin;$i++){
		my ($chr,$pos,$ori,$sline)=@{$new_sin[$i]};
		my @tmp;
		
		my $scount=0;
		foreach my $id(split(/,/,$sline)){
			$scount++;
			push @tmp,$infoarray[$id];
		}
		my $line=join(",",@tmp);
		
		my $count=scalar(split(/,/,$line));
		
		print OUT "$i\t$chr\t$pos\t$ori\t$count\t$scount\t$line\n";
	}
	close OUT;
}

if(1){
	my @infoarray;
	# my @countarray;
	open IN,"$infile" or die $!;
	while(<IN>){
		chomp;
		my @t=split;
		$infoarray[$t[0]]=[($t[4],$t[5])];
	}
	close IN;
	
	# open IN,"$infile.raw" or die $!;
	# while(<IN>){
		# chomp;
		# my @t=split;
		# $countarray[$t[0]]=$t[5];
	# }
	# close IN;
	
	open OUT,">$outfile" or die $!;
	open OUTA,">$outfile.all" or die $!;
	
	@sv=sort {
		$a->[0] <=>	 $b->[0] or $a->[1] <=> $b->[1];
	}@sv;
	
	#ida idb chra posa oria chrb posb orib count id bar jcount
	# my $R = Statistics::R->new();
	# my $cmd=qq{
# library(MASS)
# a = read.table("$gapfile")
# x<-a[,1]
# x<-x[which(x>$rlen)]
# x<-x[which(x<$mlen)]
# subx<-x[which(x>(quantile(x,$exs)))]
# res=fitdistr(subx, "lognormal")
# Ngap=length(x);
# };
	# $R->run($cmd);
	my @out;
	my $last=-1;
	foreach my $ref (@sv){
		if($ref->[0] != $last){
			if(@out > 0){
				my $line=join(",",@out);
				print OUT "$last\t$line\n";
			}
			@out=();
			$last=$ref->[0];
		}
		
		my $sup1=$infoarray[$ref->[0]][0];
		my $sup2=$infoarray[$ref->[1]][0];
		my $m1=$infoarray[$ref->[0]][1];
		my $m2=$infoarray[$ref->[1]][1];
		my $sup3=$ref->[8];
		my $m3=$ref->[11];
		
		my @q;
		# my (@id1,@id2);
		# my $c1=0;
		# my $c2=0;
		# foreach my $k(split(/,/,$ref->[9])){
			# my ($i,$j)=split(/-/,$k);
			# push @id1,$i;
			# push @id2,$j;
		# }
		
		# foreach my $i(@id1){
			# if($countarray[$i] >= $low){
				# $c1=1;
			# }
		# }
		
		# foreach my $i(@id2){
			# if($countarray[$i] >= $low){
				# $c2=1;
			# }
		# }
		# unless ($c1 and $c2){
			# push @q,"LowQual1";
		# }
		
		my $mergelow1=int(0.8*$m1*$low);
		my $mergelow2=int(0.8*$m2*$low);
		# if($sup1 <= $mergelow1 or $sup2 <= $mergelow2 or $m1 > $Nmerge*2*0.8 or $m2 > $Nmerge*2*0.8){
		if($sup1 <= $mergelow1 or $sup2 <= $mergelow2){
			push @q,"LowQual1";
		}
		
		my $mergelows=(int($m3*0.35*$lows) > $lows ? int($m3*0.35*$lows):$lows);
		if($sup3 <= $mergelows or $m3 > ((2*$Nmerge)**2) * 0.4){
			push @q,"LowQual2";
		}
		
		my $qual;
		if(@q>0){
			$qual=join(",",@q);
		}else{
			$qual="PASS";
		}
		
		
		print OUTA join("\t",@{$ref}[0,1],$qual,@{$ref}[2..11]),"\n";
		
		if($qual eq "PASS"){
			push @out,$ref->[1].":".$ref->[8];
		}
	}
	
	if(@out > 0){
		my $line=join(",",@out);
		print OUT "$last\t$line\n";
	}
	
	close OUT;
	close OUTA;
}

sub searchBest{
	my ($i,$j,$info,$share,$exist,$ref)=@_;
	
	my $last="";
	my %temp;
	my $checkid;
	$temp{"$i-$j"}=$share->[$i]->{$j}->[0];
	while(1){
		my @id=sort {$temp{$b} <=> $temp{$a} or $a cmp $b} keys %temp;
		my $name=join("-",@id);
		if ($last eq $name){
			$checkid=$id[0];
			last;
		}
		$last=$name;
		my $k=$id[0];
		my ($x,$y)=split(/-/,$k);
		%temp=();
		$temp{"$x-$y"}=$share->[$x]->{$y}->[0];
		&singleBest($info,$share,$exist,\%temp,$x,$y);
	}
	
	($i,$j)=split(/-/,$checkid);
	%temp=();
	$temp{"$i-$j"}=$share->[$i]->{$j}->[0];
	&singlelowBest($info,$share,$exist,\%temp,$i,$j);
	
	
	my ($oria,$orib);
	my ($chra,$chrb);
	my (@posa,@posb);
	my (@id,%bar);
	my $num;
	# my $max=-1;
	foreach my $key (sort {$temp{$b} <=> $temp{$a} or $a cmp $b} keys %temp){
		my ($x,$y)=split(/-/,$key);
		my ($chr1,$pos1,$ori1)=@{$info}[4*$x,4*$x+1,4*$x+2];
		my ($chr2,$pos2,$ori2)=@{$info}[4*$y,4*$y+1,4*$y+2];
		unless(defined $oria){
			$oria=$ori1;
		}
		
		unless(defined $orib){
			$orib=$ori2;
		}
		
		unless(defined $chra){
			$chra=$chr1;
		}
		
		unless(defined $chrb){
			$chrb=$chr2;
		}
		
		if($temp{$key} >= $lows){
			push @posa,$pos1;
			push @posb,$pos2;
		}
		# if($max > 0){
			# if($max == $temp{$key}){
				# push @posa,$pos1;
				# push @posb,$pos2;
			# }
		# }else{
			# push @posa,$pos1;
			# push @posb,$pos2;
			# $max=$temp{$key};
		# }

		push @id,$key;
		foreach my $b(split(/,/,$share->[$x]->{$y}->[1])){
			$bar{$b}=1;
		}
		$exist->{$x}->{$y}=1;
	}
	
	@posa=sort{$a <=> $b} @posa;
	
	@posb=sort{$a <=> $b} @posb;
	
	push @{$ref},$chra;
	if($oria eq "R"){
		push @{$ref},$posa[-1];
	}else{
		push @{$ref},$posa[0];
	}
	
	push @{$ref},$oria;
	
	push @{$ref},$chrb;
	if($orib eq "R"){
		push @{$ref},$posb[-1];
	}else{
		push @{$ref},$posb[0];
	}
	
	push @{$ref},$orib;
	
	$num= keys %bar;

	push @{$ref},$num;
	
	push @{$ref},join(",",@id);
	
	push @{$ref},join(",",(keys %bar));
	
	return;
}

sub singleBest{
	my ($info,$share,$exist,$ref,$i,$j)=@_;
	my ($chr1,$pos1,$ori1,$count1)=@{$info}[4*$i,4*$i+1,4*$i+2,4*$i+3];
	my ($chr2,$pos2,$ori2,$count2)=@{$info}[4*$j,4*$j+1,4*$j+2,4*$j+3];
	my ($start,$end);
	
	if($ori1 eq "R"){
		$start=$pos1-($Nmerge*$bin)+$bin;
		$end=$pos1+($Nmerge*$bin);
	}else{
		$start=$pos1-($Nmerge*$bin);
		$end=$pos1+($Nmerge*$bin)-$bin;
	}
	
	my ($min_i,$max_i)=($i,$i);
	while(1){
		$min_i--;
		if($min_i <0){
			$min_i++;
			last;
		}
		my ($chr,$pos,$ori)=@{$info}[4*$min_i,4*$min_i+1,4*$min_i+2];
		if($chr ne $chr1){
			$min_i++;
			last;
		}
		
		if($ori ne $ori1){
			$min_i++;
			last;
		}
		
		if($pos < $start){
			$min_i++;
			last;
		}	
	}
	
	while(1){
		$max_i++;
		if($max_i >= ($#{$info}+1)/4){
			$max_i--;
			last;
		}
		my ($chr,$pos,$ori)=@{$info}[4*$max_i,4*$max_i+1,4*$max_i+2];
		if($chr ne $chr1){
			$max_i--;
			last;
		}
		
		if($ori ne $ori1){
			$max_i--;
			last;
		}
		
		if($pos > $end){
			$max_i--;
			last;
		}	
	}
	
	
	if($ori2 eq "R"){
		$start=$pos2-($Nmerge*$bin)+$bin;
		$end=$pos2+($Nmerge*$bin);
	}else{
		$start=$pos2-($Nmerge*$bin);
		$end=$pos2+($Nmerge*$bin)-$bin;
	}
	
	my ($min_j,$max_j)=($j,$j);
	while(1){
		$min_j--;
		if($min_j <0){
			$min_j++;
			last;
		}
		my ($chr,$pos,$ori)=@{$info}[4*$min_j,4*$min_j+1,4*$min_j+2];
		if($chr ne $chr2){
			$min_j++;
			last;
		}
		
		if($ori ne $ori2){
			$min_j++;
			last;
		}
		
		if($pos < $start){
			$min_j++;
			last;
		}	
	}
	
	while(1){
		$max_j++;
		if($max_j >= ($#{$info}+1)/4){
			$max_j--;
			last;
		}
		my ($chr,$pos,$ori)=@{$info}[4*$max_j,4*$max_j+1,4*$max_j+2];
		if($chr ne $chr2){
			$max_j--;
			last;
		}
		
		if($ori ne $ori2){
			$max_j--;
			last;
		}
		
		if($pos > $end){
			$max_j--;
			last;
		}	
	}
	
	for(my $x=$min_i;$x<=$max_i;$x++){
		for(my $y=$min_j;$y<=$max_j;$y++){
			if($x< $y){
				if(exists $share->[$x]->{$y}){
					next if exists $exist->{$x}->{$y};
					next unless $info->[4*$x+3] >= $low;
					next unless $info->[4*$y+3] >= $low;
					next unless $share->[$x]->{$y}->[0] >= $lows;
					$ref->{"$x-$y"}=$share->[$x]->{$y}->[0];
				}
			}
		}
	}
	
	return;
}

sub singlelowBest{
	my ($info,$share,$exist,$ref,$i,$j)=@_;
	my ($chr1,$pos1,$ori1,$count1)=@{$info}[4*$i,4*$i+1,4*$i+2,4*$i+3];
	my ($chr2,$pos2,$ori2,$count2)=@{$info}[4*$j,4*$j+1,4*$j+2,4*$j+3];
	my ($start,$end);
	
	if($ori1 eq "R"){
		$start=$pos1-($Nmerge*$bin)+$bin;
		$end=$pos1+($Nmerge*$bin);
	}else{
		$start=$pos1-($Nmerge*$bin);
		$end=$pos1+($Nmerge*$bin)-$bin;
	}
	
	my ($min_i,$max_i)=($i,$i);
	while(1){
		$min_i--;
		if($min_i <0){
			$min_i++;
			last;
		}
		my ($chr,$pos,$ori)=@{$info}[4*$min_i,4*$min_i+1,4*$min_i+2];
		if($chr ne $chr1){
			$min_i++;
			last;
		}
		
		if($ori ne $ori1){
			$min_i++;
			last;
		}
		
		if($pos < $start){
			$min_i++;
			last;
		}	
	}
	
	while(1){
		$max_i++;
		if($max_i >= ($#{$info}+1)/4){
			$max_i--;
			last;
		}
		my ($chr,$pos,$ori)=@{$info}[4*$max_i,4*$max_i+1,4*$max_i+2];
		if($chr ne $chr1){
			$max_i--;
			last;
		}
		
		if($ori ne $ori1){
			$max_i--;
			last;
		}
		
		if($pos > $end){
			$max_i--;
			last;
		}	
	}
	
	
	if($ori2 eq "R"){
		$start=$pos2-($Nmerge*$bin)+$bin;
		$end=$pos2+($Nmerge*$bin);
	}else{
		$start=$pos2-($Nmerge*$bin);
		$end=$pos2+($Nmerge*$bin)-$bin;
	}
	
	my ($min_j,$max_j)=($j,$j);
	while(1){
		$min_j--;
		if($min_j <0){
			$min_j++;
			last;
		}
		my ($chr,$pos,$ori)=@{$info}[4*$min_j,4*$min_j+1,4*$min_j+2];
		if($chr ne $chr2){
			$min_j++;
			last;
		}
		
		if($ori ne $ori2){
			$min_j++;
			last;
		}
		
		if($pos < $start){
			$min_j++;
			last;
		}	
	}
	
	while(1){
		$max_j++;
		if($max_j >= ($#{$info}+1)/4){
			$max_j--;
			last;
		}
		my ($chr,$pos,$ori)=@{$info}[4*$max_j,4*$max_j+1,4*$max_j+2];
		if($chr ne $chr2){
			$max_j--;
			last;
		}
		
		if($ori ne $ori2){
			$max_j--;
			last;
		}
		
		if($pos > $end){
			$max_j--;
			last;
		}	
	}
	
	for(my $x=$min_i;$x<=$max_i;$x++){
		for(my $y=$min_j;$y<=$max_j;$y++){
			if($x< $y){
				if(exists $share->[$x]->{$y}){
					next if exists $exist->{$x}->{$y};
					$ref->{"$x-$y"}=$share->[$x]->{$y}->[0];
				}
			}
		}
	}
	
	return;
}






