use strict;
use warnings;

die "Usage:$0 <cluster file> <seg link> <id link> <split link file>\n" if @ARGV !=4;

my ($cluster,$seg,$id,$out)=@ARGV;

my %barset;
open IN,"$cluster" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	my $id=$t[0];
	my @data;
	&getbar($t[-1],\@data);
	$barset{$id}=[@data];
}

close IN;

my %segset;
open IN,"$seg" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	my @link=split(/,/,$t[-1]);
	if(@link== 2){
		$link[0]=~ /(\d+)-/;
		my $id1=$1;
		$link[1]=~ /(\d+)-/;
		my $id2=$1;
		
		if(exists $segset{$id1}){
			my $check=0;
			foreach my $curid(@{$segset{$id1}}){
				$check =1 if $curid == $id2;
			}
			push @{$segset{$id1}},$id2 if $check == 0;
			
		}else{
			push @{$segset{$id1}},$id2;
		}
		
		if(exists $segset{$id2}){
			my $check=0;
			foreach my $curid(@{$segset{$id2}}){
				$check =1 if $curid == $id1;
			}
			push @{$segset{$id2}},$id1 if $check == 0;
			
		}else{
			push @{$segset{$id2}},$id1;
		}
	}
}
close IN;


open IN,"$id" or die $!;
open OUT,">$out" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	my $baseid=$t[0];
	my @segid;
	if(exists $segset{$baseid}){
		@segid=@{$segset{$baseid}};
	}
	
	my @caseid;
	my @possseg;
	foreach my $n(split(/,/,$t[1])){
		$n=~ /(\S+):/;
		my $curid=$1;
		if(@segid >0){
			my $check=0;
			foreach my $m (@segid){
				$check =1 if $m == $curid;
			}
			if ($check == 1){
				push @possseg,$curid;
				next;
			}
		}
		push @caseid,$curid;
	}
	
	my $segline;
	if(@possseg >0){
		$segline= join(",",@possseg);
	}else{
		$segline = "NULL";
	}
	
	my @overlap;
	while(@caseid > 0){
		if(@overlap == 0){
			@overlap=&getover($barset{$baseid},$barset{$caseid[0]});
		}else{
			my @cur;
			for (my $i=0;$i<=$#caseid;$i++){
				#get init share
				my @share=&getover($barset{$baseid},$barset{$caseid[$i]});
				#compare share-share
				my @tmp=&getover([@overlap],[@share]);
				
				my $len1=@share;
				my $len2=@tmp;
				
				if(($len2/$len1) > 0.75 ){
					my $line="$caseid[$i]:$len1:";
					$line.=sprintf("%.4f",($len2/$len1));
					push @cur,$line;
					
					splice(@caseid,$i,1);
				}	
				
			}
			print OUT "$baseid\t$segline\t",join(",",@cur),"\n";
			@overlap=();
		}
	}
}
close IN;
close OUT;


sub getbar{
	my $line=$_[0];
	my $ref=$_[1];
	my @t=split(/,/,$line);
	foreach my $b(@t){
		$b=~ /^(\d+)-/;
		push @{$ref},$1;
	}
	return;
}


sub getover{
	my ($base,$case)=@_;
	my @bar;
	my %tmphash;
	foreach my $n(@{$base}){
		$tmphash{$n}=1;
	}
	
	foreach my $n (@{$case}){
		push @bar,$n if exists $tmphash{$n};
	}
	
	return @bar;
}
