use strict;
use warnings;

die "Usage: $0 <region mark file> <split id file> <link chain file> <merge size> <merge limit>\n" unless @ARGV==5;

my ($filter,$split,$out,$mergesize,$mergemax)=@ARGV;

my %seginfo;
my %sup;

open IN,"$split" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	my @sid=split(/,/,$t[1]);
	foreach my $i(@sid){
		$seginfo{$t[0]}{$i}=1;
	}
	
	my @lid=split(/,/,$t[2]);
	
	foreach my $i(@lid){
		my @nt=split(/:/,$i);
		$sup{$t[0]}{$nt[0]}=$nt[1];
		$sup{$nt[0]}{$t[0]}=$nt[1];
	}
}
close IN;

my %info;
my %pos;
my %checkid;
my %clusterinfo;
my %svinfo;
my @svid;
my %chaintoid;
my %idtochain;
my $chainid=0;
my $sid=0;

open IN,"$filter" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	$info{$t[0]}=[@t];
	push @{$pos{$t[3]}},[($t[4],$t[0])];
	push @{$pos{$t[5]}},[($t[6],$t[0])];
	
	next unless $t[10] eq "PASS";
	$clusterinfo{$t[1]}{$t[0]}=1;
	$clusterinfo{$t[2]}{$t[0]}=1;
	$svinfo{$t[0]}=[($t[1],$t[2])];
	push @svid,$t[0];
	$checkid{$t[0]}=1;	
}
close IN;

foreach my $i(@svid){
	next if exists $idtochain{$i};
	my @chainset;
	&searchchain($i,\@chainset);
	foreach my $line(@chainset){
		$chainid++;
		$chaintoid{$chainid}=[@{$line}];
		foreach my $j(@{$line}){
			push @{$idtochain{$j}},$chainid;
		}
	}
}

open OUT,">$out" or die $!;
print OUT "EventID\tSvID\tBreakID1\tBreakID2\tChrA\tPosA\tChrB\tPosB\tShareBarcode\tRealType\tSimpleType\tQualityScore\tComprehensiveFilter\tHeatmapFilter\tPhaseFilter\tMapQFilter\tPairEndFilter\tLocalizationFilter\tBlackList\tControlList\tSegmentCheck1\tSegmentCheck2\tSVchain\n";
while(scalar(keys %checkid) > 0){
	my @tmp= sort chrtopos (keys %checkid);
	my $case=shift @tmp;
	my @merge;
	&mergeSV($case,\@merge);
	@merge=sort chrtopos @merge;
	
	$sid++;
	foreach my $id (@merge){
		my @data;
		if(exists $checkid{$id}){
			delete $checkid{$id};
			@data=@{$info{$id}};
			splice(@data,7,0,$sup{$data[1]}{$data[2]});
			unshift @data,"S$sid";
		}else{
			@data=@{$info{$id}};
			splice(@data,7,0,$sup{$data[1]}{$data[2]});
			unshift @data,"*S$sid";
		}
		
		my %chains;
		foreach my $chainid(@{$idtochain{$id}}){
			my $ref=$chaintoid{$chainid};
			my $line=join("-",@{$ref});
			$chains{$line}=1;
		}
		
		my $oneline=join(",",(keys %chains));
		if(keys %chains == 0){
			push @data,"NULL";
		}else{
			push @data,$oneline;
		}
		
		if($data[16]=~ /PASS/){
			$data[16]=~ /:(\d+)-(\d+)/;
			$data[5]=$1;
			$data[7]=$2;
		}
		
		unless(@merge > $mergemax){
			print OUT join("\t",@data),"\n";
		}
	}
}
close OUT;


%info=();
%pos=();
%checkid=();
%clusterinfo=();
%svinfo=();
@svid=();
%chaintoid=();
%idtochain=();
$chainid=0;
$sid=0;

open IN,"$filter" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	$info{$t[0]}=[@t];
	push @{$pos{$t[3]}},[($t[4],$t[0])];
	push @{$pos{$t[5]}},[($t[6],$t[0])];
	
	next unless $t[10]=~ /PASS/;
	# next if $t[9]=~ /BAD_REGION/;
	$clusterinfo{$t[1]}{$t[0]}=1;
	$clusterinfo{$t[2]}{$t[0]}=1;
	$svinfo{$t[0]}=[($t[1],$t[2])];
	push @svid,$t[0];
	$checkid{$t[0]}=1;	
}
close IN;

foreach my $i(@svid){
	next if exists $idtochain{$i};
	my @chainset;
	&searchchain($i,\@chainset);
	foreach my $line(@chainset){
		$chainid++;
		$chaintoid{$chainid}=[@{$line}];
		foreach my $j(@{$line}){
			push @{$idtochain{$j}},$chainid;
		}
	}
}


open OUT,">$out.NoRegionFilter" or die $!;
print OUT "EventID\tSvID\tBreakID1\tBreakID2\tChrA\tPosA\tChrB\tPosB\tShareBarcode\tRealType\tSimpleType\tQualityScore\tComprehensiveFilter\tHeatmapFilter\tPhaseFilter\tMapQFilter\tPairEndFilter\tLocalizationFilter\tBlackList\tControlList\tSegmentCheck1\tSegmentCheck2\tSVchain\n";
while(scalar(keys %checkid) > 0){
	my @tmp= sort chrtopos (keys %checkid);
	my $case=shift @tmp;
	my @merge;
	&mergeSV($case,\@merge);
	@merge=sort chrtopos @merge;
	
	$sid++;
	foreach my $id (@merge){
		my @data;
		if(exists $checkid{$id}){
			delete $checkid{$id};
			@data=@{$info{$id}};
			splice(@data,7,0,$sup{$data[1]}{$data[2]});
			unshift @data,"S$sid";
		}else{
			@data=@{$info{$id}};
			splice(@data,7,0,$sup{$data[1]}{$data[2]});
			unshift @data,"*S$sid";
		}
		
		my %chains;
		foreach my $chainid(@{$idtochain{$id}}){
			my $ref=$chaintoid{$chainid};
			my $line=join("-",@{$ref});
			$chains{$line}=1;
		}
		
		my $oneline=join(",",(keys %chains));
		if(keys %chains == 0){
			push @data,"NULL";
		}else{
			push @data,$oneline;
		}
		
		if($data[16]=~ /PASS/){
			$data[16]=~ /:(\d+)-(\d+)/;
			$data[5]=$1;
			$data[7]=$2;
		}
		
		unless(@merge > $mergemax){
			print OUT join("\t",@data),"\n";
		}
	}
}
close OUT;

sub mergeSV{
	my $id=$_[0];
	my $ref=$_[1];
	my %all;
	$all{$id}=1;
	
	my @data=@{$info{$id}};
	my ($chr1,$pos1)=($data[3],$data[4]);
	my ($chr2,$pos2)=($data[5],$data[6]);
	foreach my $ref (@{$pos{$chr1}}){
		my ($ppos,$pid)=@{$ref};
		if(abs($ppos-$pos1) <= $mergesize){
			next if exists $all{$pid};
			next if $pid == $id;
			$all{$pid}=1;
		}
	}
			
	foreach my $ref (@{$pos{$chr2}}){
		my ($ppos,$pid)=@{$ref};
		if(abs($ppos-$pos2) <= $mergesize){
			next if exists $all{$pid};
			next if $pid == $id;
			$all{$pid}=1;
		}
	}
	
	@{$ref}=keys %all;
	return;
}

sub searchchain{
	my $baseid=$_[0];
	my $result=$_[1];
	my @conchain;
	
	my (@tmp1,@tmp2);
	my (%hash1,%hash2);
	my ($b,$f)=@{$svinfo{$baseid}};
	if(exists $seginfo{$f}){
		foreach my $i(keys %{$seginfo{$f}}){
			if(exists $clusterinfo{$i}){
				foreach my $j(keys %{$clusterinfo{$i}}){
					$hash1{$j}=1;
				}
			}
		}
	}
	@tmp1=keys %hash1;
	
	if(exists $seginfo{$b}){
		foreach my $i(keys %{$seginfo{$b}}){
			if(exists $clusterinfo{$i}){
				foreach my $j(keys %{$clusterinfo{$i}}){
					$hash2{$j}=1;
				}
			}
		}
	}
	@tmp2=keys %hash2;
	
	if(@tmp1 >0 and @tmp2 > 0){
		foreach my $i(@tmp1){
			foreach my $j(@tmp2){
				push @conchain,[($j,$baseid,$i)];
			}
		}
	}elsif(@tmp1 > 0){
		foreach my $i(@tmp1){
			push @conchain,[($baseid,$i)];
		}
	}elsif(@tmp2 >0){
		foreach my $i(@tmp2){
			push @conchain,[($i,$baseid)];
		}
	}else{
		push @conchain,[($baseid)];
		@{$result}=@conchain;
		return;
		
	}
	
	my $edge=1;
	my $n=0;
	while($edge){
		$n++;
		my $end_count=0;
		my @tmp_conchain;
		
		foreach my $ref(@conchain){
			@tmp1=();
			@tmp2=();
			%hash1=();
			%hash2=();
			my %existsnode;
			foreach my $node(@{$ref}){
				$existsnode{$node}=1;
			}
			
			#forword
			my $checkid=$ref->[-1];
			my ($b,$f)=@{$svinfo{$checkid}};
			
			if(exists $seginfo{$f}){
				foreach my $i(keys %{$seginfo{$f}}){
					if(exists $clusterinfo{$i}){
						foreach my $j(keys %{$clusterinfo{$i}}){
							next if exists $existsnode{$j};
							$hash1{$j}=1;
						}
					}
				}
			}
			
			if(exists $seginfo{$b}){
				foreach my $i(keys %{$seginfo{$b}}){
					if(exists $clusterinfo{$i}){
						foreach my $j(keys %{$clusterinfo{$i}}){
							next if exists $existsnode{$j};
							$hash1{$j}=1;
						}
					}
				}
			}
			@tmp1=keys %hash1;
			#backword
			$checkid=$ref->[0];
			($b,$f)=@{$svinfo{$checkid}};
			
			if(exists $seginfo{$f}){
				foreach my $i(keys %{$seginfo{$f}}){
					if(exists $clusterinfo{$i}){
						foreach my $j(keys %{$clusterinfo{$i}}){
							next if exists $existsnode{$j};
							$hash2{$j}=1;
						}
					}
				}
			}
			
			if(exists $seginfo{$b}){
				foreach my $i(keys %{$seginfo{$b}}){
					if(exists $clusterinfo{$i}){
						foreach my $j(keys %{$clusterinfo{$i}}){
							next if exists $existsnode{$j};
							$hash2{$j}=1;
						}
					}
				}
			}
			@tmp2=keys %hash2;
			
			if(@tmp1 >0 and @tmp2 > 0){
				foreach my $i(@tmp1){
					foreach my $j(@tmp2){
						my @tmp_line=@{$ref};
						push @tmp_line,$i;
						unshift @tmp_line,$j;
						push @tmp_conchain,[@tmp_line];
					}
				}
			}elsif(@tmp1 > 0){
				foreach my $i(@tmp1){
					my @tmp_line=@{$ref};
					push @tmp_line,$i;
					push @tmp_conchain,[@tmp_line];
				}
			}elsif(@tmp2 >0){
				foreach my $i(@tmp2){
					my @tmp_line=@{$ref};
					unshift @tmp_line,$i;
					push @tmp_conchain,[@tmp_line];
				}
			}else{
				my @tmp_line=@{$ref};
				push @tmp_conchain,[@tmp_line];
				$end_count++;
			}	
		}

		if($end_count == scalar @conchain){
			$edge=0;
		}else{
			@conchain=@tmp_conchain;
		}
		
	}
	
	@{$result}=@conchain;
	return;
}

sub chrtopos{
	my @ta=@{$info{$a}};
	my @tb=@{$info{$b}};
	my $nchra1=$ta[3];
	my $nchra2=$ta[5];
	$nchra1=~ s/chr//;
	$nchra2=~ s/chr//;
	$nchra1 =23 if $nchra1 eq "X";
	$nchra1 =24 if $nchra1 eq "Y";
	$nchra1 =25 if $nchra1 eq "M";
	$nchra2 =23 if $nchra2 eq "X";
	$nchra2 =24 if $nchra2 eq "Y";
	$nchra2 =25 if $nchra2 eq "M";
	my $posa1=$ta[4];
	my $posa2=$ta[6];
	
	my $nchrb1=$tb[3];
	my $nchrb2=$tb[5];
	$nchrb1=~ s/chr//;
	$nchrb2=~ s/chr//;
	$nchrb1 =23 if $nchrb1 eq "X";
	$nchrb1 =24 if $nchrb1 eq "Y";
	$nchrb1 =25 if $nchrb1 eq "M";
	$nchrb2 =23 if $nchrb2 eq "X";
	$nchrb2 =24 if $nchrb2 eq "Y";
	$nchrb2 =25 if $nchrb2 eq "M";
	my $posb1=$tb[4];
	my $posb2=$tb[6];
	
	if($nchra1=~ /^\d+$/ and $nchrb1=~ /^\d+$/ ){
		if($nchra2=~ /^\d+$/ and $nchrb2=~ /^\d+$/ ){
			$nchra1 <=> $nchrb1 or $posa1 <=> $posb1 or $nchra2 <=> $nchrb2 or $posa2 <=> $posb2;
		}else{
			$nchra1 <=> $nchrb1 or $posa1 <=> $posb1 or $nchra2 cmp $nchrb2 or $posa2 <=> $posb2;
		}
	}else{
		if($nchra2=~ /^\d+$/ and $nchrb2=~ /^\d+$/ ){
			$nchra1 cmp $nchrb1 or $posa1 <=> $posb1 or $nchra2 <=> $nchrb2 or $posa2 <=> $posb2;
		}else{
			$nchra1 cmp $nchrb1 or $posa1 <=> $posb1 or $nchra2 cmp $nchrb2 or $posa2 <=> $posb2;
		}
	}
}
