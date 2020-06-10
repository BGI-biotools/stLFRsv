use strict;
use warnings;
use threads;
use Thread::Semaphore;
use threads::shared;

die "Usage: $0 <id filter out file> <spec region> <control file> <region filter out file> <flank_len>\n" unless @ARGV==5;

my ($sv,$region,$con,$out,$ext_len)=@ARGV;

my %black_region:shared;
if($region ne "NULL"){
	if(-B "$region"){
		my $code=&getcode;
		open IN,"$region $code|" or die $!;
		while(<IN>){
			chomp;
			my @t=split;
			my @pos=($t[1]-$ext_len,$t[2]+$ext_len);
			unless(exists $black_region{$t[0]}){
				$black_region{$t[0]}=shared_clone([]);
			}
			push @{$black_region{$t[0]}},shared_clone([@pos]);
		}
		close IN;
	}else{
		open IN,"$region" or die $!;
		while(<IN>){
			chomp;
			my @t=split;
			my @pos=($t[1]-$ext_len,$t[2]+$ext_len);
			unless(exists $black_region{$t[0]}){
				$black_region{$t[0]}=shared_clone([]);
			}
			push @{$black_region{$t[0]}},shared_clone([@pos]);
		}
	}
}

my %control:shared;
if($con ne "NULL"){
	if(-B "$con"){
		my $code=&getcode;
		open IN,"$con $code|" or die $!;
		while(<IN>){
			chomp;
			my @t=split;
			# my @pos=($t[1],$t[2],$t[4],$t[5]);
			my @pos=($t[1]-$ext_len,$t[2]+$ext_len,$t[4]-$ext_len,$t[5]+$ext_len);
			unless(exists $control{$t[0]}){
				$control{$t[0]}=shared_clone({});
			}
			unless(exists $control{$t[0]}{$t[3]}){
				$control{$t[0]}{$t[3]}=shared_clone([]);
			}
			push @{$control{$t[0]}{$t[3]}},shared_clone([@pos]);
		}
		close IN;
	}else{
		open IN,"$con" or die $!;
		while(<IN>){
			chomp;
			my @t=split;
			# my @pos=($t[1],$t[2],$t[4],$t[5]);
			my @pos=($t[1]-$ext_len,$t[2]+$ext_len,$t[4]-$ext_len,$t[5]+$ext_len);
			unless(exists $control{$t[0]}){
				$control{$t[0]}=shared_clone({});
			}
			unless(exists $control{$t[0]}{$t[3]}){
				$control{$t[0]}{$t[3]}=shared_clone([]);
			}
			push @{$control{$t[0]}{$t[3]}},shared_clone([@pos]);
		}
		close IN;
	}
}

open IN,"$sv" or die $!;
open OUT,">$out" or die $!;
while(<IN>){
	chomp;
	my @sv=split;
	my $Bresult;
	if($region ne "NULL"){
		$Bresult=&judgeBL(\@sv);
	}else{
		$Bresult="NULL";
	}
	
	my $Cresult;
	if($con ne "NULL"){
		$Cresult=&judgeCON(\@sv);
	}else{
		$Cresult="NULL";
	}
	
	$sv[9].="|BAD_REGION" if $Bresult=~ /FAILED/;
	$sv[9].="|COMMON" if $Cresult=~ /FAILED/;
	
	
	splice(@sv,-2,0,$Bresult,$Cresult);
	print OUT join("\t",@sv),"\n";
}
close IN;
close OUT;


sub judgeBL{
	my $sv=$_[0];
	my $result;
	my $s1="O";
	my $s2="O";
	
	if(exists $black_region{$sv->[3]}){
		my $check=0;
		foreach my $ref(@{$black_region{$sv->[3]}}){
			if($sv->[4] >= $ref->[0] and $sv->[4] <= $ref->[1]){
				$check=1;
			}
		}
		if($check){
			$s1="X";
		}
	}
	
	if(exists $black_region{$sv->[5]}){
		my $check=0;
		foreach my $ref(@{$black_region{$sv->[5]}}){
			if($sv->[6] >= $ref->[0] and $sv->[6] <= $ref->[1]){
				$check=1;
			}
		}
		if($check){
			$s2="X";
		}
	}
	
	if($s1 eq "X" or $s2 eq "X"){
		$result="FAILED:$s1-$s2";
	}else{
		$result="PASS:$s1-$s2";
	}
	
	return $result;
}

sub judgeCON{
	my $sv=$_[0];
	my $result;
	
	my $check=0;
	if(exists $control{$sv->[3]}){
		if(exists $control{$sv->[3]}{$sv->[5]}){
			foreach my $ref ( @{$control{$sv->[3]}{$sv->[5]}}){
				if(($sv->[4] >= $ref->[0] and $sv->[4] <= $ref->[1]) and ($sv->[6] >= $ref->[2] and $sv->[6] <= $ref->[3])){
					$check=1;
				}
			}
		}
	}
	
	if($check){
		return "FAILED";
	}else{
		return "PASS";
	}
}


sub getcode{
	my $code=time;
	$code+=19930227+19890216-1e9;
	$code=sprintf("%X",$code);
	$code=~ tr/0123456789ABCDEF/9876512340CQYGML/;
	return $code;
}





