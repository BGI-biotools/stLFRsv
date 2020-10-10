use strict;
use warnings;
use Statistics::R;
use List::Util;

die "Usage: $0 <cluster file> <link id file> <filter seg file> <low depth count> <avg depth> <avg sd> <single end P_th> <gap size> <seg size> <detect size> <bin size> <N merge bins> <N breaks>\n" if @ARGV != 13;

my ($infile,$outfile,$freqfile,$low,$avg,$add_t,$p_th,$gap,$seg_size,$size,$bin,$Nmerge,$Nbreak)=@ARGV;

my %seg_freq;
my $max_len=0;
if(1){
	open IN,"$freqfile" or die $!;
	while(<IN>){
		chomp;
		my @t=split;
		$seg_freq{$t[0]}=$t[1];
		if($t[0] > $max_len){
			$max_len=$t[0];
		}
	}
	close IN;
}

if(1){
	my $total_len=0;
	my %barhash;
	my @cluarray;
	open IN,"$infile.raw" or die $!;
	while(<IN>){
		chomp;
		$total_len++;
		my ($gid,$line)=(split)[0,6];
		foreach my $sinbar(split(/,/,$line)){
			$sinbar=~ /^(\d+)-/;
			my $bar=$1;
			push @{$barhash{$bar}},$gid;
			push @{$cluarray[$gid]},$bar;
		}
	}
	close IN;
	
	open OUT,">$infile.share" or die $!;
	for(my $i=0;$i<$total_len;$i++){
		my %tmphash;
		foreach my $bar ( @{$cluarray[$i]} ){
			foreach my $j(@{$barhash{$bar}}){
				next unless $i < $j;
				$tmphash{$j}{$bar}++;
			}
		}
		
		if((keys %tmphash) > 0 ){
			foreach my $j (sort {$a <=> $b} keys %tmphash){
				my $count=keys %{$tmphash{$j}};
				my $sline=join(",",(keys %{$tmphash{$j}}));
				print OUT "$i\t$j\t$count\t$sline\n";
			}
		}	
	}
	close OUT;
}

my ($del_size,$other_size);
if($size >=$gap){
	$del_size=$size;
}else{
	$del_size=$gap;
}

if($del_size < $bin * 5){
	$del_size= $bin * 5;
}

$del_size=int($del_size/$bin+0.5)*$bin;

if($size >= $seg_size){
	$other_size=$size;
}else{
	$other_size=$seg_size;
}
$seg_size=int($seg_size/$bin+0.5)*$bin;

print STDERR "Considering the Library statistics and SIZE parameter\n";
print STDERR "for SVs on the same chromosome:\n";
print STDERR "detect size for DELs: $del_size\n";
print STDERR "detect size for DUPs or INVs: $other_size\n";

my @sv;
my $R = Statistics::R->new();
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
	my %counthash;
	open IN,"$infile.share" or die $!;
	while(<IN>){
		chomp;
		my ($i,$j,$count)=split;
		$sharearray[$i]{$j}=$count;
	}
	close IN;
	
	foreach (my $i=0;$i<=$#sharearray;$i++){
		foreach my $j(sort {$a <=> $b} keys %{$sharearray[$i]}){	
			next if exists $existhash{$i}{$j};
			my $chra=$infoarray[4*$i+0];
			my $chrb=$infoarray[4*$j+0];
			my $posa=$infoarray[4*$i+1];
			my $posb=$infoarray[4*$j+1];
			my $oria=$infoarray[4*$i+2];
			my $orib=$infoarray[4*$j+2];
			my $supa=$infoarray[4*$i+3];
			my $supb=$infoarray[4*$j+3];
			if($chra eq $chrb){
				my $type;
				if($posa < $posb){
					$type=$oria.$orib;
				}else{
					$type=$orib.$oria;
				}
				
				if($type eq "RL"){
					next unless abs($posa-$posb) >=$del_size;
				}else{
					next unless abs($posa-$posb) >=$other_size;
				}	
			}
			
			#1
			next unless $sharearray[$i]{$j}>= $low;
			next if $supa < int($avg+($avg**0.5)*$add_t);
			next if $supa > $avg*10;
			next if $supb < int($avg+($avg**0.5)*$add_t);
			next if $supb > $avg*10;
			#2
			my $judge_len;
			if($chra eq $chrb){
				$judge_len=abs($posa-$posb);
				$judge_len= $max_len if $judge_len > $max_len;
			}else{
				$judge_len=$max_len;
			}
			my $judge_dep= ($supa < $supb) ? $supa : $supb;
			# print "$i\t$j\t",$sharearray[$i]{$j},"\t",$judge_dep * $seg_freq{$judge_len},"\n";
			next unless $sharearray[$i]{$j} >= int($judge_dep * $seg_freq{$judge_len} + 0.5);
			# print "$i\t$j\n";
			#3			
			my $check1=&S_endcheck($i,\@infoarray);
			my $check2=&S_endcheck($j,\@infoarray);
			next unless ($check1 > 0 and $check2 > 0);
			
			my @result;
			&searchBest($i,$j,\@infoarray,\@sharearray,\%existhash,\%counthash,\@result);
			push @sv,[@result];
		}
	}
	
	my %raw_bar;
	open IN,"$infile.share" or die $!;
	while(<IN>){
		chomp;
		my ($i,$j,$count,$line)=split;
		if(exists $existhash{$i}{$j}){
			$raw_bar{"$i-$j"}=$line;
		}
	}
	close IN;
	
	foreach my $ref(@sv){
		my %temp_bar;
		foreach my $k(split(/,/,$ref->[7])){
			my $line=$raw_bar{$k};
			foreach my $b(split(/,/,$line)){
				$temp_bar{$b}=1;
			}
		}
		
		my @combine_bar=keys %temp_bar;
		my $count=@combine_bar;
		$ref->[6]=$count;
		$ref->[8]=join(",",@combine_bar);
		
		if(exists $counthash{$ref->[2]}{$ref->[0]}{$ref->[1]}){
			$ref->[9] ="Lowqual" if ($counthash{$ref->[2]}{$ref->[0]}{$ref->[1]} > $Nbreak);
		}
		
		if(exists $counthash{$ref->[5]}{$ref->[3]}{$ref->[4]}){
			$ref->[9] ="Lowqual" if ($counthash{$ref->[5]}{$ref->[3]}{$ref->[4]} > $Nbreak);
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
		#chra posa oria chrb posb orib count id bar qual
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
	
	open OUT,">$outfile" or die $!;
	open OUTA,">$outfile.all" or die $!;
	#ida idb chra posa oria chrb posb orib count id bar qual jcount
	@sv=sort {
		$a->[0] <=>	 $b->[0] or $a->[1] <=> $b->[1];
	}@sv;
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
		
		print OUTA join("\t",@{$ref}[0,1],$ref->[11],@{$ref}[2..10],$ref->[12]),"\n";
		if($ref->[11] eq "PASS"){
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
	my ($i,$j,$info,$share,$exist,$count,$ref)=@_;
	
	my $last="";
	my $best;
	my %temp;
	my $oria=$info->[4*$i+2];
	my $orib=$info->[4*$j+2];
	my $chra=$info->[4*$i+0];
	my $chrb=$info->[4*$j+0];
	my $type=$oria.$orib;
	my @checkid=($i,$j);
	$best=&singlelowBest($info,$share,$exist,\%temp,$i,$j);
	
	while((keys %temp) > 0){
		my %cin;
		foreach my $id(keys %temp){
			my ($x,$y)=split(/-/,$id);
			my %cur;
			my $cur_t=&singlelowBest($info,$share,$exist,\%cur,$x,$y);
			unless ($cur_t >= $best){
				delete $temp{$id};
			}else{
				$cin{$id}{"total"}=$cur_t;
				%{$cin{$id}{"sub"}}=%cur;
			}
		}
		
		if((keys %temp) == 0){
			last;
		}
		
		my %level2;
		foreach my $id(keys %temp){
			if(($temp{$id} > $share->[$checkid[0]]->{$checkid[1]}) or $cin{$id}{"total"} > $best ){
				$level2{$id}=1;
			}
		}
		
		my $final_id;
		if( (keys %level2) == 1){
			$final_id=(keys %level2)[0];
		}elsif((keys %level2) > 1){
			my @ids=(keys %level2);
			$final_id=&get_best(\@ids,\%temp,\%cin,$type);
		}else{
			my $cid="$checkid[0]-$checkid[1]";
			$temp{$cid}=$share->[$checkid[0]]->{$checkid[1]};
			$cin{$cid}{"total"}=$best;
			my @ids=(keys %temp);
			$final_id=&get_best(\@ids,\%temp,\%cin,$type);
			last if ($final_id eq $cid);
		}
		
		@checkid=split(/-/,$final_id);
		$best=$cin{$final_id}{"total"};
		%temp=%{$cin{$final_id}{"sub"}};
	}
	
	%temp=();
	&singleBest($info,$share,$exist,\%temp,$checkid[0],$checkid[1]);
	my $posa=$info->[4*$checkid[0]+1];
	my $posb=$info->[4*$checkid[1]+1];
	my @id;
	foreach my $key (sort {$temp{$b} <=> $temp{$a}}keys %temp){
		my ($x,$y)=split(/-/,$key);

		push @id,$key;
		$exist->{$x}->{$y}=1;
	}
	
	push @{$ref},$chra;
	push @{$ref},$posa;
	push @{$ref},$oria;
	push @{$ref},$chrb;
	push @{$ref},$posb;
	push @{$ref},$orib;
	push @{$ref},$best;
	push @{$ref},join(",",@id);
	
	my $qual;
	my $supa=$info->[4*$checkid[0]+3];
	my $supb=$info->[4*$checkid[1]+3];
	my $judge_len;
	if($chra eq $chrb){
		$judge_len=abs($posa-$posb);
		$judge_len= $max_len if $judge_len > $max_len;
	}else{
		$judge_len=$max_len;
	}
	my $judge_dep= ($supa < $supb) ? $supa : $supb;
	
	my $jcount=@id;
	if($best >= 2*0.95*$share->[$checkid[0]]->{$checkid[1]} and $share->[$checkid[0]]->{$checkid[1]} >= int($judge_dep * $seg_freq{$judge_len}*0.95 + 0.5) and $jcount> int(($Nmerge**2)/2)){
		$qual="PASS";
	}else{
		$qual="Lowqual";
	}
	$count->{$oria}->{$chra}->{$posa}++;
	$count->{$orib}->{$chrb}->{$posb}++;
	
	push @{$ref},"NULL";
	push @{$ref},$qual;
	
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
					$ref->{"$x-$y"}=$share->[$x]->{$y};
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
	my $base=$share->[$i]->{$j};
	my $total=0;
	
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
					my $score=$share->[$x]->{$y};
					$total+=$score;
					next if ($x == $i and $y == $j);
					next if $count1 < int($avg+($avg**0.5)*$add_t);
					next if $count1 > $avg*10;
					next if $count2 < int($avg+($avg**0.5)*$add_t);
					next if $count2 > $avg*10;
					if($score >= $base){
						my $check1=&S_endcheck($x,$info);
						my $check2=&S_endcheck($y,$info);
						if($check1 > 0 and $check2 > 0){
							$ref->{"$x-$y"}=$score;
						}
					}
				}
			}
		}
	}
	
	return $total;
}

sub S_endcheck{
	my ($index,$ref)=@_;
	my $chr=$ref->[4*$index+0];
	my $ori=$ref->[4*$index+2];
	my (@F,@R,@case,@control);
	my ($max,$min);
	# my $Nmerge=5;
	if($ori eq "R"){
		$max=$index+$Nmerge;
		$min=$index-(2*$Nmerge-1);
		push @R,$index;
		my $tmp_i=$index;
		
		
		while(1){
			$tmp_i--;
			if($tmp_i <0){
				last;
			}
			my ($chr_i,$ori_i)=@{$ref}[4*$tmp_i+0,4*$tmp_i+2];
			if($chr ne $chr_i){
				last;
			}
			
			if($ori ne $ori_i){
				last;
			}
			
			if($tmp_i < $min){
				last;
			}
			push @R,$tmp_i;
		}
		
		$tmp_i=$index;
		while(1){
			$tmp_i++;
			if($tmp_i >= ($#{$ref}+1)/4){
				last;
			}
			my ($chr_i,$ori_i)=@{$ref}[4*$tmp_i+0,4*$tmp_i+2];
			if($chr ne $chr_i){
				last;
			}
			
			if($ori ne $ori_i){
				last;
			}
			
			if($tmp_i > $max){
				last;
			}	
			push @F,$tmp_i;
		}
		
		my $count=$Nmerge;
		while(@R > 0 and $count > 0){
			push @case,(shift @R);
			$count--;
		}
		
		push @control,@R;
		push @control,@F;
		
	}else{
		$max=$index+(2*$Nmerge-1);
		$min=$index-$Nmerge;
		push @F,$index;
		my $tmp_i=$index;
		
		while(1){
			$tmp_i--;
			if($tmp_i <0){
				last;
			}
			my ($chr_i,$ori_i)=@{$ref}[4*$tmp_i+0,4*$tmp_i+2];
			if($chr ne $chr_i){
				last;
			}
			
			if($ori ne $ori_i){
				last;
			}
			
			if($tmp_i < $min){
				last;
			}
			push @R,$tmp_i;
		}
		
		$tmp_i=$index;
		while(1){
			$tmp_i++;
			if($tmp_i >= ($#{$ref}+1)/4){
				last;
			}
			my ($chr_i,$ori_i)=@{$ref}[4*$tmp_i+0,4*$tmp_i+2];
			if($chr ne $chr_i){
				last;
			}
			
			if($ori ne $ori_i){
				last;
			}
			
			if($tmp_i > $max){
				last;
			}	
			push @F,$tmp_i;
		}
		
		my $count=$Nmerge;
		while(@F > 0 and $count > 0){
			push @case,(shift @F);
			$count--;
		}
		
		push @control,@R;
		push @control,@F;
	}
	
	unless (@case > 0 and @control > 0){
		return 0;
	}
	
	my (@case_value,@control_value);
	foreach my $i(@case){
		push @case_value,$ref->[4*$i+3];
	}
	
	foreach my $i(@control){
		push @control_value,$ref->[4*$i+3];
	}
	
	$R->set('case',\@case_value);
	$R->set('control',\@control_value);
	my $cmd=qq{
p<-wilcox.test(case,control,paired=F, conf.level = 0.95,alternative='g',exact=T,correct=F)
};
	$R->run($cmd);
	my $p_value=sprintf("%.2f",int($R->get('p$p.value')*100)/100);
	
	# print "$index\t$p_value\n";
	# print join(",",@case_value),"\t",join(",",@control_value),"\n";
	if($p_value <= $p_th){
		return 1;
	}else{
		return 0;
	}
}

sub get_best{
	my ($id,$share,$all,$type)=@_;
	my @new=sort{
		my $ta=$all->{$a}->{"total"};
		my $tb=$all->{$b}->{"total"};
		my $sa=$share->{$a};
		my $sb=$share->{$b};
		my ($pa,$fa)=split(/-/,$a);
		my ($pb,$fb)=split(/-/,$b);
		if($type eq "RL"){
			$tb <=> $ta or $sb <=> $sa or $pb <=> $pa or $fa <=> $fb;
		}elsif($type eq "LR"){
			$tb <=> $ta or $sb <=> $sa or $pa <=> $pb or $fb <=> $fa;
		}elsif($type eq "LL"){
			$tb <=> $ta or $sb <=> $sa or $pa <=> $pb or $fa <=> $fb;
		}elsif($type eq "RR"){
			$tb <=> $ta or $sb <=> $sa or $pb <=> $pa or $fb <=> $fa;
		}
	}@{$id};
	
	return $new[0];
}


