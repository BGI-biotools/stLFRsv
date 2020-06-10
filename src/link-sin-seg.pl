use strict;
use warnings;

die "Usage: $0 <sbf file> <cluster file> <seg link file>\n" if @ARGV != 3;

my ($seg_file,$clu_file,$out_file)=@ARGV;

my %barhash;
my %cluhash;
open IN,"$clu_file" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	my $id=$t[0];
	my $ori=$t[3];
	my @data;
	&getbar($t[-1],\@data);
	
	foreach my $c (@data){
		my $bar=$c->[0];
		my $index=$c->[1];
		push @{$barhash{$bar}{$index}},[($id,$ori)];
	}
}

close IN;

my $seg_index=$seg_file.".bfi";
open INDEX,"$seg_index" or die $!;
binmode(INDEX);
my $buf;
my %index;
while(read(INDEX,$buf,8)){
	my $bar=unpack("%64Q",$buf);
	read(INDEX,$buf,8);
	my $off=unpack("%64Q",$buf);
	$index{$bar}=$off;
}
close INDEX;

open IN,"$seg_file" or die $!;
binmode(IN);
open OUT,">$out_file" or die $!;
foreach my $bar (keys %barhash){
	my $temp=$bar >> 20;
	seek(IN,$index{$temp},0);
	while(read(IN,$buf,8)){
		my $curbar=unpack("%64Q",$buf);
		read(IN,$buf,4);
		my $index=unpack("%32L",$buf);
	
		read(IN,$buf,32);
		my $name=unpack("Z32",$buf);
	
		read(IN,$buf,4);
		my $s=unpack("%32L",$buf);
	
		read(IN,$buf,4);
		my $e=unpack("%32L",$buf);
	
		read(IN,$buf,4);
		my $sup=unpack("%32L",$buf);
		
		last if ($curbar >> 20) != $temp;
		last if ($curbar & 0xFFFFF) > ($bar & 0xFFFFF);
		next if ($curbar & 0xFFFFF) < ($bar & 0xFFFFF);
		if(exists $barhash{$bar}{$index}){
			my @tmp;
			foreach my $ref (@{$barhash{$bar}{$index}}){
				push @tmp,$ref->[0]."-".$ref->[1];	
			}
			print OUT "$bar\t$index\t$name\t$s\t$e\t$sup\t".join(",",@tmp)."\n";
		}else{
			print OUT "$bar\t$index\t$name\t$s\t$e\t$sup\tNULL\n";
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
		$b=~ /^(\d+)-(\d+):/;
		push @{$ref},[($1,$2)];
	}
	return;
}
