use strict;
use warnings;

die "Usage: $0 <judge sv file> <split id file> <filter out file> <sameid count>\n" unless @ARGV==4;

my ($judge,$split,$out,$sth)=@ARGV;

my %info;
my %id1;
my %id2;
my %result1;
my %result2;

open IN,"$judge" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	$info{$t[0]}=[@t];
	if($t[9] eq "PASS"){
		push @{$id1{$t[1]}},$t[0];
		push @{$id1{$t[2]}},$t[0];
	
		$id2{"$t[1]-$t[2]"}=$t[0];
		$id2{"$t[2]-$t[1]"}=$t[0];
	}
}
close IN;

foreach my $key(keys %id1){
	if(@{$id1{$key}} > $sth){
		foreach my $id(@{$id1{$key}}){
			$result1{$id}=1;
		}
	}
}

open IN,"$split" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	if($t[-1]=~ /,/){
		my %id_set;
		my @data=split(/,/,$t[-1]);
		foreach my $i(@data){
			$i=~ /^(\d+):(\d+):/;
			my $cid=$1;
			my $num=$2;
			if(exists $id2{"$t[0]-$cid"}){
				$id_set{$id2{"$t[0]-$cid"}}=$num;
			}
		}
		
		if((keys %id_set) >1){
			my $max=0;
			my $mid="";
			foreach my $i(keys %id_set){
				if($id_set{$i} > $max){
					$max=$id_set{$i};
					$mid=$i;
				}
			}
			
			foreach my $i(keys %id_set){
				next if $i == $mid;
				$result2{$i}=1;
			}
		}
	}
}
close IN;

open OUT,">$out" or die $!;
foreach my $i(sort {$a <=> $b} keys %info){
	my @data=@{$info{$i}};
	my $check=0;
	
	if(exists $result1{$i}){
		if($data[9] eq "PASS"){
			$check=1;
			push @data,"FAILED";
		}else{
			push @data,"SKIPED";
		}
	}else{
		if($data[9] eq "PASS"){
			push @data,"PASS";
		}else{
			push @data,"SKIPED";
		}
		
	}
	
	if(exists $result2{$i}){
		if($data[9] eq "PASS"){
			$check=1;
			push @data,"FAILED";
		}else{
			push @data,"SKIPED";
		}
	}else{
		if($data[9] eq "PASS"){
			push @data,"PASS";
		}else{
			push @data,"SKIPED";
		}
	}
	
	if($check == 1 and $data[9] eq "PASS"){
		$data[9] = "FAILED";
	}
	
	print OUT join("\t",@data,"\n");
}



