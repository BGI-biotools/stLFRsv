use strict;
use warnings;
use Bio::DB::HTS;
use Statistics::Descriptive;
use threads;
use Thread::Semaphore;
use threads::shared;

die "Usage: $0 <cluster file> <split id file> <indexed bam file> <phase info dir> <out file> <check size> <ncpu>\n" if @ARGV != 7;

# my ($cluster,$split,$bamfile,$phase_dir,$region,$con,$outfile,$size,$ncpu)=@ARGV;
my ($cluster,$split,$bamfile,$phase_dir,$outfile,$size,$ncpu)=@ARGV;

my %clusterinfo:shared;
open IN,"$cluster" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	$clusterinfo{$t[0]}=shared_clone([@t]);
}
close IN;

my %svinfo:shared;
my %svresult:shared;
my %link;
my $svid=0;
open IN,"$split" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	my $id0=$t[0];
	my @toid=split(/,/,$t[2]);
	&get_type($id0,\@toid);
}
close IN;

my %phaseinfo:shared;
my %phaseregion:shared;
if($phase_dir ne "NULL"){
	opendir DIR,"$phase_dir" or die $!;
	foreach my $i(readdir DIR){
		next unless $i=~ /barcode\.phase$/;
		open IN,"$phase_dir/$i" or die $!;
		my $name=$i;
		$name=~ s/\.barcode\.phase//;
		while(<IN>){
			chomp;
			my @t=split;
			unless(exists $phaseinfo{$name}){
				$phaseinfo{$name}=shared_clone({});
			}
			$phaseinfo{$name}{$t[0]}=$t[-1];
		}
		close IN;
		
		open IN,"$phase_dir/$name.region" or die $!;
		my $c=0;
		while(<IN>){
			chomp;
			my @t=split;
			$c++;
			unless(exists $phaseregion{$name}){
				$phaseregion{$name}=shared_clone({});
			}
			$phaseregion{$name}{$c}=shared_clone([@t]);
		}
		close IN;
	}
	closedir DIR;
}

# my %black_region:shared;
# if($region ne "NULL"){
	# if(-B "$region"){
		# my $code=&getcode;
		# open IN,"$region $code|" or die $!;
		# while(<IN>){
			# chomp;
			# my @t=split;
			# my @pos=($t[1],$t[2]);
			# unless(exists $black_region{$t[0]}){
				# $black_region{$t[0]}=shared_clone([]);
			# }
			# push @{$black_region{$t[0]}},shared_clone([@pos]);
		# }
		# close IN;
	# }else{
		# open IN,"$region" or die $!;
		# while(<IN>){
			# chomp;
			# my @t=split;
			# my @pos=($t[1],$t[2]);
			# unless(exists $black_region{$t[0]}){
				# $black_region{$t[0]}=shared_clone([]);
			# }
			# push @{$black_region{$t[0]}},shared_clone([@pos]);
		# }
	# }
# }

# my %control:shared;
# if($con ne "NULL"){
	# if(-B "$con"){
		# my $code=&getcode;
		# open IN,"$con $code|" or die $!;
		# while(<IN>){
			# chomp;
			# my @t=split;
			# my @pos=($t[1],$t[2],$t[4],$t[5]);
			# unless(exists $control{$t[0]}){
				# $control{$t[0]}=shared_clone({});
			# }
			# unless(exists $control{$t[0]}{$t[3]}){
				# $control{$t[0]}{$t[3]}=shared_clone([]);
			# }
			# push @{$control{$t[0]}{$t[3]}},shared_clone([@pos]);
		# }
		# close IN;
	# }else{
		# open IN,"$con" or die $!;
		# while(<IN>){
			# chomp;
			# my @t=split;
			# my @pos=($t[1],$t[2],$t[4],$t[5]);
			# unless(exists $control{$t[0]}){
				# $control{$t[0]}=shared_clone({});
			# }
			# unless(exists $control{$t[0]}{$t[3]}){
				# $control{$t[0]}{$t[3]}=shared_clone([]);
			# }
			# push @{$control{$t[0]}{$t[3]}},shared_clone([@pos]);
		# }
		# close IN;
	# }
# }

my $semaphore=new Thread::Semaphore($ncpu);
open OUT,">$outfile" or die $!;
foreach my $svid(sort {$a <=> $b}keys %svinfo){
	$semaphore->down();
	my $thread=threads->new(\&processOne,$svid);
	foreach my $t (threads->list(threads::joinable)){
		my $value=$t->join();
		if($value){
			exit(1);
		}
	}
}

&waitquit;

foreach my $svid(sort {$a <=> $b}keys %svresult){
	print OUT $svresult{$svid}."\n";
}
close OUT;

# sub testF{
	# sleep(10);
	# $semaphore->up();
	# return;
# }

sub processOne{
	my $svid=$_[0];
	my $sv=$svinfo{$svid};
	my $typeone;
	if($sv->[6] eq "LL" or $sv->[6] eq "RR"){
		$typeone="F";
	}elsif($sv->[6] eq "LR" or $sv->[6] eq "RL"){
		$typeone="R";
	}
	
	my @result;
	&checkfour($sv,$typeone,\@result);
	
	# my ($P,$F)=(0,0);
	# foreach my $i(@result){
		# $P++ if $i=~ /PASS/;
		# $F++ if $i=~ /FAILED/;
	# }
	
	# if($region eq "NULL"){
		# $P--;
	# }
	
	# if($con eq "NULL"){
		# $P--;
	# }
	
	my $final;
	# if($F >=2){
		# $final="FAILED";
	# }elsif($F == 1 and $P >= 2){
		# $final="PASS";
	# }elsif($F ==0 and $P >=1){
		# $final="PASS";
	# }else{
		# $final="FAILED";
	# }
	
	if($result[0]=~ /PASS/ and $result[1]=~ /PASS/){
		$final="PASS";
	}elsif($result[0]=~ /SUSPECT/ and $result[1]=~ /PASS/){
		$final="PASS";
	}elsif($result[0]!~ /FAILED/ and $result[1]!~ /FAILED/ and $result[2]=~ /PASS/){
		$final="PASS";
	}elsif($result[1]=~ /PASS/ and $result[2]=~ /PASS/){
		$final="PASS";
	}elsif($result[0]=~ /PASS/ and $result[2]=~ /PASS/){
		$final="PASS";
	}else{
		$final="FAILED";
	}
	
	my $simple_type;
	if($sv->[2] eq $sv->[4]){
		if($sv->[6] eq "RL"){
			$simple_type="DEL";
		}elsif($sv->[6] eq "LR"){
			$simple_type="DUP";
		}elsif($sv->[6] eq "LL"){
			$simple_type="INV1";
		}elsif($sv->[6] eq "RR"){
			$simple_type="INV2";
		}
	}else{
		if($sv->[6] eq "RL"){
			$simple_type="TRA1";
		}elsif($sv->[6] eq "LR"){
			$simple_type="TRA2";
		}elsif($sv->[6] eq "LL"){
			$simple_type="TRA3";
		}elsif($sv->[6] eq "RR"){
			$simple_type="TRA4";
		}
	}
	my $line=join("\t",$svid,@{$sv}[0..6],$simple_type,$final,@result);
	
	$svresult{$svid}=$line;
	$semaphore->up();
	return;
}

sub checkfour{
	my ($ref,$type,$out)=@_;
	
	##1 heatmap
	my (@bar1,@bar2);
	if($type eq "F"){
		@bar1=&share_dis($ref->[2],$ref->[3],$ref->[4],$ref->[5],"R",10,10,$size);
		@bar2=&share_dis($ref->[4],$ref->[5],$ref->[2],$ref->[3],"R",10,10,$size);
	}else{
		@bar1=&share_dis($ref->[2],$ref->[3],$ref->[4],$ref->[5],"L",10,10,$size);
		@bar2=&share_dis($ref->[4],$ref->[5],$ref->[2],$ref->[3],"R",10,10,$size);
	}
	
	my @sb;
	$type=$ref->[6];
	foreach my $i(0..$#bar1){
		my $n;
		if ($bar1[$i] eq "N" or $bar2[$i] eq "N"){
			$n="N";
		}else{
			$n=$bar1[$i]*0.707106781+$bar2[$i]*0.707106781;
		}
		push @sb,$n;
	}
	my $peak=&searchbin(\@sb);
	
	my (@binR,@binL,@binS1,@binS2);
	my (@barR,@barL,@barS1,@barS2);
	
	push @binR,$peak + 3;
	push @binR,$peak + 4;
	push @binR,$peak + 5;
	
	push @binL,$peak - 3;
	push @binL,$peak - 4;
	push @binL,$peak - 5;
	
	push @binS1,$peak -1;
	push @binS1,$peak -2;
	push @binS1,$peak -3;
	push @binS1,$peak -4;
	push @binS1,$peak -5;

	
	push @binS2,$peak +1;
	push @binS2,$peak +2;
	push @binS2,$peak +3;
	push @binS2,$peak +4;
	push @binS2,$peak +5;
	
	my $judge1="NULL";
	foreach my $i(@binR){
		goto PASS if $sb[$i] eq "N";
		push @barR,$sb[$i];
	}
	
	foreach my $i(@binL){
		goto PASS if $sb[$i] eq "N";
		push @barL,$sb[$i];
	}
	
	foreach my $i(@binS1){
		goto PASS if $sb[$i] eq "N";
		push @barS1,$sb[$i];
	}
	
	foreach my $i(@binS2){
		goto PASS if $sb[$i] eq "N";
		push @barS2,$sb[$i];
	}
	
	if($sb[$peak] < (0.3* $ref->[7]) ){
		$judge1="LOW";
		goto PASS;
	}
	
	my $stat = Statistics::Descriptive::Full->new();
	
	#1
	$stat->clear();
	$stat->add_data(@barR);
	my $s1=$stat->median();
	my $delta1=$sb[$peak]-$s1;
	my $ratio1=$delta1/$sb[$peak];
	
	#2
	$stat->clear();
	$stat->add_data(@barL);
	my $s2=$stat->median();
	my $delta2=$sb[$peak]-$s2;
	my $ratio2=$delta2/$sb[$peak];
	
	#3
	my $sum;
	for my $i (0..$#barS1){
		$sum+=$barS1[$i]-$barS2[$i];
	}
	$sum=abs($sum)/($#barS1+1);
	
	my $delta3=$sum;

	my $ratio3=$sum/$sb[$peak];
	
	my ($j1,$j2,$j3);
	
	if($ratio1 <= 0.25){
		$j1=1;
	}elsif($ratio1 >0.25 and $ratio1 <=0.5){
		$j1=2;
	}elsif($ratio1 > 0.5){
		$j1=3;
	}
	
	if($ratio2 <= 0.25){
		$j2=1;
	}elsif($ratio2 >0.25 and $ratio2 <=0.5){
		$j2=2;
	}elsif($ratio2 > 0.5){
		$j2=3;
	}
	
	if($ratio3 <= 0.1){
		$j3=1;
	}elsif($ratio3 >0.1 and $ratio3 <=0.2){
		$j3=2;
	}elsif($ratio3 > 0.2){
		$j3=3;
	}
	
	if($ref->[2] eq $ref->[4]){
		if($type eq "RL"){
			if($j1 ==3 and $j2 <=2 and $j3 ==3){
				$judge1="PASS";
			}elsif($j1 ==3 and $j2 ==1 and $j3 >=2){
				$judge1="PASS";
			}elsif($j1 ==3 and $j3 >=2){
				$judge1="SUSPECT";
			}else{
				$judge1="FAILED";
			}
		}
		
		if($type eq "LR"){
			if($j1 <=2 and $j2 ==3 and $j3 ==3){
				$judge1="PASS";
			}elsif($j1 ==1 and $j2 ==3 and $j3 >=2){
				$judge1="PASS";
			}elsif($j2 ==3 and $j3 >=2){
				$judge1="SUSPECT";
			}else{
				$judge1="FAILED";
			}
		}
		
		if($type eq "RR"){
			if($j1 ==3 and $j2 <=2 and $j3 ==3){
				$judge1="PASS";
			}elsif($j1 ==3 and $j2 ==1 and $j3 >=2){
				$judge1="PASS";
			}elsif(abs($j1-$j2)<=1 and $j1 >=2 and $j2 >=2 and $j3 == 1 ){
				$judge1="PASS";
			}elsif(abs($j1-$j2)<=1 and $j1 >=2 and $j2 >=2 and $ratio3 <= 0.5 ){
				$judge1="SUSPECT";
			}else{
				$judge1="FAILED";
			}
		}
		
		if($type eq "LL"){
			if($j1 <=2 and $j2 ==3 and $j3 ==3){
				$judge1="PASS";
			}elsif($j1 ==1 and $j2 ==3 and $j3 >=2){
				$judge1="PASS";
			}elsif(abs($j1-$j2)<=1 and $j1 >=2 and $j2 >=2 and $j3 == 1 ){
				$judge1="PASS";
			}elsif(abs($j1-$j2)<=1 and $j1 >=2 and $j2 >=2 and $ratio3 <= 0.5){
				$judge1="SUSPECT";
			}else{
				$judge1="FAILED";
			}
		}
	
	}else{
		if($type eq "RL" or $type eq "RR"){
			if($j1 ==3 and $j2 <=2 and $j3 ==3){
				$judge1="PASS";
			}elsif($j1 ==3 and $j2 ==1 and $j3 >=2){
				$judge1="PASS";
			}elsif(abs($j1-$j2)<=1 and $j1 >=2 and $j2 >=2 and $j3 == 1 ){
				$judge1="PASS";
			}elsif(abs($j1-$j2)<=1 and $j1 >=2 and $j2 >=2 and $ratio3 <= 0.5 ){
				$judge1="SUSPECT";
			}else{
				$judge1="FAILED";
			}
		}
		
		if($type eq "LR" or $type eq "LL"){
			if($j1 <=2 and $j2 ==3 and $j3 ==3){
				$judge1="PASS";
			}elsif($j1 ==1 and $j2 ==3 and $j3 >=2){
				$judge1="PASS";
			}elsif(abs($j1-$j2)<=1 and $j1 >=2 and $j2 >=2 and $j3 == 1 ){
				$judge1="PASS";
			}elsif(abs($j1-$j2)<=1 and $j1 >=2 and $j2 >=2 and $ratio3 <= 0.5 ){
				$judge1="SUSPECT";
			}else{
				$judge1="FAILED";
			}
		}
	}

	PASS:
	my $line1;
	if($judge1 eq "LOW"){
		$line1="FAILED:NA(NA)|NA(NA)|NA(NA)";
	}elsif($judge1 ne "NULL"){
		$delta1=sprintf("%.4f",$delta1);
		$delta2=sprintf("%.4f",$delta2);
		$delta3=sprintf("%.4f",$delta3);
		$ratio1=sprintf("%.4f",$ratio1);
		$ratio2=sprintf("%.4f",$ratio2);
		$ratio3=sprintf("%.4f",$ratio3);
		
		$line1="$judge1:$ratio1($delta1)|$ratio2($delta2)|$ratio3($delta3)";
		
	}else{
		$line1="$judge1:NA(NA)|NA(NA)|NA(NA)";
	}
	
	$out->[0]=$line1;
	

	#barcode phase HP
	my $judge2=&judgeshare($ref);
	$out->[1]=$judge2;
	
	#mapQ
	my $judge3=&judgemap($ref);
	$out->[2]=$judge3;

	#region
	# my $judge4=&judgeBL($ref);
	# $out->[3]=$judge4;
	
	#control
	# my $judge5=&judgeCON($ref);
	# $out->[4]=$judge5;

	return;
}

# sub judgeBL{
	# my $sv=$_[0];
	# my $result;
	# my $s1="O";
	# my $s2="O";
	
	# if(exists $black_region{$sv->[2]}){
		# my $check=0;
		# foreach my $ref(@{$black_region{$sv->[2]}}){
			# if($sv->[3] >= $ref->[0] and $sv->[3] <= $ref->[1]){
				# $check=1;
			# }
		# }
		# if($check){
			# $s1="X";
		# }
	# }
	
	# if(exists $black_region{$sv->[4]}){
		# my $check=0;
		# foreach my $ref(@{$black_region{$sv->[4]}}){
			# if($sv->[5] >= $ref->[0] and $sv->[5] <= $ref->[1]){
				# $check=1;
			# }
		# }
		# if($check){
			# $s2="X";
		# }
	# }
	
	# if($s1 eq "X" or $s2 eq "X"){
		# $result="FAILED:$s1-$s2";
	# }else{
		# $result="PASS:$s1-$s2";
	# }
	
	# return $result;
# }

# sub judgeCON{
	# my $sv=$_[0];
	# my $result;
	
	# my $check=0;
	# if(exists $control{$sv->[2]}){
		# if(exists $control{$sv->[2]}{$sv->[4]}){
			# foreach my $ref ( @{$control{$sv->[2]}{$sv->[4]}}){
				# if(($sv->[3] >= $ref->[0] and $sv->[3] <= $ref->[1]) and ($sv->[5] >= $ref->[2] and $sv->[5] <= $ref->[3])){
					# $check=1;
				# }
			# }
		# }
	# }
	
	# if($check){
		# return "FAILED";
	# }else{
		# return "PASS";
	# }
# }

sub judgemap{
	my $sv=$_[0];
	my $result;
	my $start=$sv->[3]-$size;
	my $end=$sv->[3]+$size;
	$start=0 if $start <0;
	my ($a1,$h1)=&get_mapq_from_bam($sv->[2],$start,$end);
	$start=$sv->[5]-$size;
	$end=$sv->[5]+$size;
	$start=0 if $start <0;
	my ($a2,$h2)=&get_mapq_from_bam($sv->[4],$start,$end);
	
	if($a1 ==0 or $a2 ==0){
		$result = "FAILED:NA-NA";
		goto PASS;
	}
	my $r1=sprintf("%.4f",$h1/$a1);
	my $r2=sprintf("%.4f",$h2/$a2);
	
	if($r1 < 0.5 and $r2 < 0.5){
		$result = "FAILED:$r1-$r2";
	}else{
		$result = "PASS:$r1-$r2";
	}
PASS:
	return $result;
}


sub judgeshare{
	my $sv=$_[0];
	
	my $result;
	#judge the region
	my $chr1=$sv->[2];
	my $pos1=$sv->[3];
	my $check=1;
	my $hp1;
	if(exists $phaseregion{$chr1}){
		foreach my $i(keys %{$phaseregion{$chr1}}){
			if($pos1 >= $phaseregion{$chr1}{$i}->[0]  and $pos1 <= $phaseregion{$chr1}{$i}->[1]){
				$check=0;
				$hp1="$chr1:$i:".$phaseregion{$chr1}{$i}->[0]."-".$phaseregion{$chr1}{$i}->[1];
				last;
			}
		}
	}
	
	if($check){
		$result="NULL:NA:NA:NA:NA:NA:NA";
		goto PASS;
	}
	
	$check=1;
	my $chr2=$sv->[4];
	my $pos2=$sv->[5];
	my $hp2;
	if(exists $phaseregion{$chr2}){
		foreach my $i(keys %{$phaseregion{$chr2}}){
			if($pos2 >= $phaseregion{$chr2}{$i}->[0]  and $pos2 <= $phaseregion{$chr2}{$i}->[1]){
				$check=0;
				$hp2="$chr2:$i:".$phaseregion{$chr2}{$i}->[0]."-".$phaseregion{$chr2}{$i}->[1];
				last;
			}	
		}
	}
	
	if($check){
		$result="NULL:NA:NA:NA:NA:NA:NA";
		goto PASS;
	}
	
	
	my @one=split(/,/,$clusterinfo{$sv->[0]}->[-1]);
	my @two=split(/,/,$clusterinfo{$sv->[1]}->[-1]);
	my %local_hash1;
	my %local_hash2;
	foreach my $i(@one){
		$i=~ /^(\d+)\-/;
		my $bar=$1;
		$local_hash1{$bar}=1;
	}
	
	my @share;
	foreach my $i(@two){
		$i=~ /^(\d+)\-/;
		my $bar=$1;
		$local_hash2{$bar}=1;
		if(exists $local_hash1{$bar}){
			push @share,$bar;
		}
	}
	
	my @stat1=(0,0,0,0);
	my @stat2=(0,0,0,0);
	my @stat3=(0,0,0,0);
	my @stat4=(0,0,0,0);
	
	foreach my $i( keys %local_hash1){
		my $bar1=$i & 0xFFFFF;
		$i =$i >> 20;
		my $bar2=$i & 0xFFFFF;
		$i =$i >> 20;
		my $bar3=$i & 0xFFFFF;
		
		my $line=join("_",$bar3,$bar2,$bar1);
	
		
		$stat1[0]+=1;
		if(exists $phaseinfo{$chr1}){
			if(exists $phaseinfo{$chr1}{$line}){
				if($phaseinfo{$chr1}{$line} eq "PAT"){
					$stat1[1]+=1;
				}else{
					$stat1[2]+=1;
				}
			}else{
				$stat1[3]+=1;
			}
		}else{
			$stat1[3]+=1;
		}
	}
	
	my $ratio1=sprintf("%.4f",$stat1[3]/$stat1[0]);
	if($ratio1 > 0.75){
		$result="NULL:UNPHASED:$ratio1:NA:NA:$hp1:$hp2";
		goto PASS;
	}
	
	
	foreach my $i( keys %local_hash2){
		my $bar1=$i & 0xFFFFF;
		$i =$i >> 20;
		my $bar2=$i & 0xFFFFF;
		$i =$i >> 20;
		my $bar3=$i & 0xFFFFF;
		
		my $line=join("_",$bar3,$bar2,$bar1);
	
		
		$stat2[0]+=1;
		if(exists $phaseinfo{$chr2}){
			if(exists $phaseinfo{$chr2}{$line}){
				if($phaseinfo{$chr2}{$line} eq "PAT"){
					$stat2[1]+=1;
				}else{
					$stat2[2]+=1;
				}
			}else{
				$stat2[3]+=1;
			}
		}else{
			$stat2[3]+=1;
		}
	}
	
	my $ratio2=sprintf("%.4f",$stat2[3]/$stat2[0]);
	if($ratio2 > 0.75){
		$result="NULL:UNPHASED:$ratio1:$ratio2:NA:$hp1:$hp2";
		goto PASS;
	}
	
	foreach my $i(@share){
		my $bar1=$i & 0xFFFFF;
		$i =$i >> 20;
		my $bar2=$i & 0xFFFFF;
		$i =$i >> 20;
		my $bar3=$i & 0xFFFFF;
		
		my $line=join("_",$bar3,$bar2,$bar1);
	
		if($chr1 eq $chr2){
			$stat3[0]+=1;
			if(exists $phaseinfo{$chr1}){
				if(exists $phaseinfo{$chr1}{$line}){
					if($phaseinfo{$chr1}{$line} eq "PAT"){
						$stat3[1]+=1;
					}else{
						$stat3[2]+=1;
					}
				}else{
					$stat3[3]+=1;
				}
			}else{
				$stat3[3]+=1;
			}
			
			my $ratio3=sprintf("%.4f",$stat3[3]/$stat3[0]);
			if($ratio3 > 0.75){
				$result="NULL:UNPHASED:$ratio1:$ratio2:$ratio3:$hp1:$hp2";
				goto PASS;
			}			
		}else{
			$stat3[0]+=1;
			if(exists $phaseinfo{$chr1}){
				if(exists $phaseinfo{$chr1}{$line}){
					if($phaseinfo{$chr1}{$line} eq "PAT"){
						$stat3[1]+=1;
					}else{
						$stat3[2]+=1;
					}
				}else{
					$stat3[3]+=1;
				}
			}else{
				$stat3[3]+=1;
			}
			my $ratio3=sprintf("%.4f",$stat3[3]/$stat3[0]);
			if($ratio3 > 0.75){
				$result="NULL:UNPHASED:$ratio1:$ratio2:$ratio3/NA:$hp1:$hp2";
				goto PASS;
			}
			
		
			$stat4[0]+=1;
			if(exists $phaseinfo{$chr2}){
				if(exists $phaseinfo{$chr2}{$line}){
					if($phaseinfo{$chr2}{$line} eq "PAT"){
						$stat4[1]+=1;
					}else{
						$stat4[2]+=1;
					}
				}else{
					$stat4[3]+=1;
				}
			}else{
				$stat4[3]+=1;
			}
			my $ratio4=sprintf("%.4f",$stat4[3]/$stat4[0]);
			if($ratio4 > 0.75){
				$result="NULL:UNPHASED:$ratio1:$ratio2:$ratio3/$ratio4:$hp1:$hp2";
				goto PASS;
			}
		}	
	}
	
	
	
	my $type;
	if($chr1 eq $chr2){
		$ratio1=sprintf("%.4f",$stat1[1]/($stat1[1]+$stat1[2]));
		if($ratio1<=0.25){
			$type="0|1";
		}elsif($ratio1 >=0.75){
			$type="1|0";
		}else{
			$type="1|1";
		}
		
		my $local_type;
		$ratio2=sprintf("%.4f",$stat2[1]/($stat2[1]+$stat2[2]));
		
		if($ratio2<=0.25){
			$local_type="0|1";
		}elsif($ratio2 >=0.75){
			$local_type="1|0";
		}else{
			$local_type="1|1";
		}
		
		if($type ne $local_type){
			$result="FAILED:MIXED:$ratio1:$ratio2:NA:$hp1:$hp2";
			goto PASS;
		}
		
		if($type eq "1|1"){
			my $dif=abs($ratio1-$ratio2);
			if($dif > 0.15){
				$result="FAILED:MIXED:$ratio1:$ratio2:NA:$hp1:$hp2";
				goto PASS;
			}
		}
		
		my $ratio3=sprintf("%.4f",$stat3[1]/($stat3[1]+$stat3[2]));
		if($ratio3<=0.25){
			$local_type="0|1";
		}elsif($ratio3 >=0.75){
			$local_type="1|0";
		}else{
			$local_type="1|1";
		}
		
		
		if($type ne $local_type){
			$result="FAILED:MIXED:$ratio1:$ratio2:$ratio3:$hp1:$hp2";
			goto PASS;
		}
		
		$result="PASS:$type:$ratio1:$ratio2:$ratio3:$hp1:$hp2";
	}else{
		$ratio1=sprintf("%.4f",$stat1[1]/($stat1[1]+$stat1[2]));
		if($ratio1<=0.25){
			$type="0|1";
		}elsif($ratio1 >=0.75){
			$type="1|0";
		}else{
			$type="1|1";
		}
		
		my $local_type;
		$ratio2=sprintf("%.4f",$stat2[1]/($stat2[1]+$stat2[2]));
		
		if($ratio2<=0.25){
			$local_type="0|1";
		}elsif($ratio2 >=0.75){
			$local_type="1|0";
		}else{
			$local_type="1|1";
		}
		
		if($type eq "1|1"){
			if($type ne $local_type){
				$result="FAILED:MIXED:$ratio1:$ratio2:NA:$hp1:$hp2";
				goto PASS;
			}
		}else{
			unless($local_type ne "1|1"){
				$result="FAILED:MIXED:$ratio1:$ratio2:NA:$hp1:$hp2";
				goto PASS;
			}
		}
		
		my $ratio3=sprintf("%.4f",$stat3[1]/($stat3[1]+$stat3[2]));
		if($ratio3<=0.25){
			$local_type="0|1";
		}elsif($ratio3 >=0.75){
			$local_type="1|0";
		}else{
			$local_type="1|1";
		}
		
		if($type eq "1|1"){
			if($type ne $local_type){
				$result="FAILED:MIXED:$ratio1:$ratio2:$ratio3/NA:$hp1:$hp2";
				goto PASS;
			}
		}else{
			unless($local_type ne "1|1"){
				$result="FAILED:MIXED:$ratio1:$ratio2:$ratio3/NA:$hp1:$hp2";
				goto PASS;
			}
		}
		
		my $ratio4=sprintf("%.4f",$stat4[1]/($stat4[1]+$stat4[2]));
		if($ratio4<=0.25){
			$local_type="0|1";
		}elsif($ratio4 >=0.75){
			$local_type="1|0";
		}else{
			$local_type="1|1";
		}
		
		if($type eq "1|1"){
			if($type ne $local_type){
				$result="FAILED:MIXED:$ratio1:$ratio2:$ratio3/$ratio4:$hp1:$hp2";
				goto PASS;
			}
		}else{
			unless($local_type ne "1|1"){
				$result="FAILED:MIXED:$ratio1:$ratio2:$ratio3/$ratio4:$hp1:$hp2";
				goto PASS;
			}
		}
		$result="PASS:$type:$ratio1:$ratio2:$ratio3/$ratio4:$hp1:$hp2";
	}
	
	PASS:
	return $result;
}



sub searchbin{
	my $ref=$_[0];
	my @id=(8,9,10,11,12);
	
	my $max=0;
	my $peak=10;
	foreach my $i(@id){
		next if $ref->[$i] eq "N";
		if( $ref->[$i] >$max){
			$peak=$i;
			$max=$ref->[$i];
		}
	}
	return $peak;
}

sub share_dis{
	my ($chr,$pos,$chrext,$posext,$type,$Lbin,$Hbin,$size)=@_;
	
	#get pos barcode
	my %checkhash;
	my $start=$pos-$size;
	my $end=$pos+$size;
	$start=0 if $start < 0;
	&get_bar_from_bam($chr,$start,$end,\%checkhash);
	
	my @region;
	if($type eq "L"){
		$start=$posext+$Lbin*$size;
		$end=$posext-$Hbin*$size;
		$end =0 if $end < 0;
		for(my $i=$start;$i>=$end;$i-=$size){
			my %tmphash;
			next if $i == 0;
			&get_bar_from_bam($chrext,$i-$size,$i,\%tmphash);
			push @region,[keys %tmphash];
		}
	}else{
		$start=$posext-$Lbin*$size;
		$end=$posext+$Hbin*$size;
		$start=0 if $start < 0;
		
		for(my $i=$start;$i<=$end;$i+=$size){
			my %tmphash;
			&get_bar_from_bam($chrext,$i,$i+$size,\%tmphash);
			push @region,[keys %tmphash];
		}
	}
	
	my @overlap;
	foreach my $ref (@region){
		my $n=0;
		foreach my $barcode (@{$ref}){
			$n++ if exists $checkhash{$barcode};
		}
		push @overlap,$n;
	}
	
	while (@overlap != ($Lbin + $Hbin +1)){
		if($type eq "L"){
			push @overlap,"N";
		}else{
			unshift @overlap,"N";
		}
	}
	

	return @overlap;
}

sub get_bar_from_bam{
	my ($chr,$start,$end,$ref)=@_;
	my $hfile = Bio::DB::HTSfile->open($bamfile);
	my $index = Bio::DB::HTSfile->index_load($hfile);
	my $header = $hfile->header_read;
	my $callback = sub {
		my ($alignment,$data) = @_;
		# my $flag=$a->flag;
		#return if $flag & 0x400;
		my $pos=($alignment->pos)+1;
		return if $pos < $start;
		return if $pos > $end;
		my $name=$alignment->qname;
		my $barcode=(split(/#/,$name))[-1];
		return if $barcode eq "0_0_0";
		$ref->{$barcode}=1;
	};
	$index->fetch($hfile,$header->parse_region("$chr:$start-$end"),$callback);
	
	return;
}

sub get_mapq_from_bam{
	my ($chr,$start,$end)=@_;
	my ($all,$high)=(0,0);
	my $hfile = Bio::DB::HTSfile->open($bamfile);
	my $index = Bio::DB::HTSfile->index_load($hfile);
	my $header = $hfile->header_read;
	my $callback = sub {
		my ($alignment,$data) = @_;
		my $flag=$alignment->flag;
		return if $flag & 0x400;
		my $pos=($alignment->pos)+1;
		return if $pos < $start;
		return if $pos > $end;
		my $q=$alignment->qual;
		$all++;
		if($q > 10){
			$high++;
		}
	};
	$index->fetch($hfile,$header->parse_region("$chr:$start-$end"),$callback);
	
	return ($all,$high);
}

sub get_type{
	my ($id,$ref)=@_;
	my @info0=@{$clusterinfo{$id}};
	foreach my $i(@{$ref}){
		my @tmp=split(/:/,$i);
		
		#jump the exists
		next if exists $link{$id}{$tmp[0]};
		
		#assign sv id
		$svid++;
		#cache the id to sv
		$link{$id}{$tmp[0]}=$svid;
		$link{$tmp[0]}{$id}=$svid;
		
		my @info1=@{$clusterinfo{$tmp[0]}};
		##creat the SV
		my @local_sv;
		if($info0[1] eq $info1[1]){
			if($info0[2] < $info1[2]){
				my $chr1=$info0[1];
				my $pos1=$info0[2];
				my $chr2=$info1[1];
				my $pos2=$info1[2];
				my $type=$info0[3].$info1[3];
				@local_sv=($id,$tmp[0],$chr1,$pos1,$chr2,$pos2,$type,$tmp[1]);
			}else{
				my $chr1=$info1[1];
				my $pos1=$info1[2];
				my $chr2=$info0[1];
				my $pos2=$info0[2];
				my $type=$info1[3].$info0[3];
				@local_sv=($tmp[0],$id,$chr1,$pos1,$chr2,$pos2,$type,$tmp[1]);
			}
		}else{
			my $nchr1=$info0[1];
			my $nchr2=$info1[1];
			$nchr1=~ s/chr//;
			$nchr2=~ s/chr//;
			$nchr1 =23 if $nchr1 eq "X";
			$nchr1 =24 if $nchr1 eq "Y";
			$nchr1 =25 if $nchr1 eq "M";
			$nchr2 =23 if $nchr2 eq "X";
			$nchr2 =24 if $nchr2 eq "Y";
			$nchr2 =25 if $nchr2 eq "M";
			
			my $check=0;
			if($nchr1=~ /^\d+$/ and $nchr2=~ /^\d+$/ ){
				if($nchr1 < $nchr2){
					$check=1;
				}
			}else{
				if($nchr1 lt $nchr2){
					$check=1;
				}
			}
			
			if($check){
				my $chr1=$info0[1];
				my $pos1=$info0[2];
				my $chr2=$info1[1];
				my $pos2=$info1[2];
				my $type=$info0[3].$info1[3];
				@local_sv=($id,$tmp[0],$chr1,$pos1,$chr2,$pos2,$type,$tmp[1]);
			
			}else{
				my $chr1=$info1[1];
				my $pos1=$info1[2];
				my $chr2=$info0[1];
				my $pos2=$info0[2];
				my $type=$info1[3].$info0[3];
				@local_sv=($tmp[0],$id,$chr1,$pos1,$chr2,$pos2,$type,$tmp[1]);
			
			}	
		}	
		$svinfo{$svid}=shared_clone([@local_sv]);
	}
	return;
}


sub getcode{
	my $code=time;
	$code+=19930227+19890216-1e9;
	$code=sprintf("%X",$code);
	$code=~ tr/0123456789ABCDEF/9876512340CQYGML/;
	return $code;
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
###############################################################






