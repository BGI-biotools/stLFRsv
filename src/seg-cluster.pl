use strict;
use warnings;
use Bio::DB::HTS;
use Statistics::R;
use Statistics::Descriptive;

die "Usage:$0 <split link file> <sin cluster file> <sbf file> <out seg ins file> <indexed bam> <filter threshold> <filter percentage> <binsize> <filter region>\n" if @ARGV !=9;

my ($split,$cluster,$sbf,$out,$bam)=@ARGV;

my $sam = Bio::DB::HTS->new(-bam  => $bam,);

my @para=split(/,/,$ARGV[5]);
my $per=$ARGV[6];
my $binsize=$ARGV[7];
my $region=$ARGV[8];

my %region1;
my $code=&getcode;
open IN,"$region $code|" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	my @pos=($t[1],$t[2]);
	push @{$region1{$t[0]}},[@pos];
}
close IN;

my $barth;
if(@para ==1){
	$barth=$para[0]*1.3;
}elsif(@para == 2){
	my $R = Statistics::R->new();
	my $cmd=qq{
library(MASS)
th<-qlnorm($per,$para[0],$para[1],lower.tail = TRUE)
};
	$R->run($cmd);
	$barth=$R->get('th');
	$R->stop();
}else{
	die "Wrong filter threshold parameter!\n";
}

#exclude the split id
my %exbar;
open IN,"$split" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	$exbar{$t[0]}=1;
}
close IN;


my %indexhash;
open IN,"$sbf.bfi" or die $!;
binmode(IN);
my $buf;
while(read(IN,$buf,8)){
	my $bar=unpack("%64Q",$buf);
	my $bar1=$bar & 0xFFFFF;
	$bar =$bar >> 20;
	my $bar2=$bar & 0xFFFFF;
	$bar =$bar >> 20;
	my $bar3=$bar & 0xFFFFF;
	read(IN,$buf,8);
	my $off=unpack("%64Q",$buf);
	
	$indexhash{$bar2}{$bar1}=$off;
}
close IN;


my $th;
open SBF,"$sbf" or die $!;
binmode(SBF);
open IN,"$cluster" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	next if exists $exbar{$t[0]};
	my $flag=0;
	my $chr=$t[1];
	my $pos=$t[2];
	$chr="chr".$chr unless ($chr=~ /chr/); 
	foreach my $ref (@{$region1{$chr}}){
		my $check1=$ref->[0];
		my $check2=$ref->[1];
		if($pos >= ($check1-$binsize)  and $pos <= ($check2+ $binsize) ){
			$flag=1;
			last;
		}
	}
	
	next if $flag == 1;
	
	my @barset;
	my %seghash;
	&getbar($t[-1],\@barset,\%seghash);
	
	$th=@barset;
	$th=int($th*0.6);
	
	my @bat;
	foreach my $bar(@barset){
		&getseg($bar,\@bat,\%seghash);
	}
	
	my @region;
	&getregion(\@region,\@bat);
	
	my @share;
	&getshare(\@region,\@barset,\@share);
	
	if (@share > 0){
		my @line;
		foreach my $ref(@share){
			push @line,$ref->[0]."-".$ref->[1]."-".$ref->[2];
		}
		print "$t[0]\t$t[1]\t$t[2]\t$t[3]\t",join(":",@line),"\n";
	}
	
	
}
close IN;




sub getbar{
	my $line=$_[0];
	my $ref=$_[1];
	my $seg=$_[2];
	my @t=split(/,/,$line);
	foreach my $b(@t){
		$b=~ /^(\d+)-/;
		push @{$ref},$1;
		$b=~ /^(\S+):/;
		$seg->{$1}=1;
	}
	return;
}

sub getseg{
	my $bar=$_[0];
	my $ref=$_[1];
	my $seg=$_[2];
	
	my $bar1=$bar & 0xFFFFF;
	$bar =$bar >> 20;
	my $bar2=$bar & 0xFFFFF;
	$bar =$bar >> 20;
	my $bar3=$bar & 0xFFFFF;
	
	my $off=$indexhash{$bar3}{$bar2};
	
	my $buf;
	seek(SBF,$off,0);
	while(read(SBF,$buf,8)){
		my $newbar=unpack("%64Q",$buf);
		read(SBF,$buf,4);
		my $index=unpack("%32L",$buf);

		read(SBF,$buf,32);
		my $name=unpack("Z32",$buf);

		read(SBF,$buf,4);
		my $s=unpack("%32L",$buf);

		read(SBF,$buf,4);
		my $e=unpack("%32L",$buf);

		read(SBF,$buf,4);
		my $sup=unpack("%32L",$buf);
		
		my $line=$newbar."-".$index;
		
		my $newbar1=$newbar & 0xFFFFF;
		$newbar =$newbar >> 20;
		my $newbar2=$newbar & 0xFFFFF;
		$newbar =$newbar >> 20;
		my $newbar3=$newbar& 0xFFFFF;
		
		next if $newbar1 < $bar1;
		next if exists $seg->{$line};
		last if $newbar1 > $bar1;
		
		push @{$ref},[($name,$s,$e)];
	}
	
	return;
}

sub getregion{
	my ($region,$bat)=@_;
	my %temphash;
	foreach my $seg(@{$bat}){
		for(my $i=int($seg->[1]/$binsize);$i<=int($seg->[2]/$binsize);$i++){
			$temphash{$seg->[0]}{$i}++;
		}
	}
	
	foreach my $chr( 
	sort {
		my $namea=$a;
		my $nameb=$b;
		$namea=~ s/^chr//;
		$nameb=~ s/^chr//;
		$namea= 23 if $namea eq "X";
		$namea= 24 if $namea eq "Y";
		$namea= 25 if $namea eq "M";
		$nameb= 23 if $nameb eq "X";
		$nameb= 24 if $nameb eq "Y";
		$nameb= 25 if $nameb eq "M";
		$namea <=> $nameb;
	}keys %temphash){
		foreach my $bin(sort {$a <=> $b} keys %{$temphash{$chr}}){
			if($temphash{$chr}{$bin} > $th){
				push @{$region},[($chr,$bin)];
			}
		}
	}
	return;
}

sub getshare{
	my ($region,$bar,$share)=@_;
	
	my %base;
	foreach my $b(@{$bar}){
		my $bar1=$b & 0xFFFFF;
		$b =$b >> 20;
		my $bar2=$b & 0xFFFFF;
		$b =$b >> 20;
		my $bar3=$b & 0xFFFFF;
		
		$base{"$bar3\_$bar2\_$bar1"}=1;
	}

	foreach my $r(@{$region}){
		my %case;
		&get_bar_from_bam($r->[0],$r->[1]*$binsize,($r->[1]+1)*$binsize,\%case);
		
		my $n=0;
		foreach my $b(keys %case){
			$n++ if exists $base{$b};
		}
		
		push @{$share},[($r->[0],$r->[1]*$binsize,$n)] if $n >= $barth;
	}
	return;
}

sub get_bar_from_bam{
	my ($chr,$start,$end,$ref)=@_;
	my @alignments = $sam->get_features_by_location(-seq_id => $chr,-start  => $start,-end => $end);
	foreach my $a(@alignments){
		my $flag=$a->flag;
		next if $flag & 0x400;
		my $pos=$a->pos;
		next if $pos < $start;
		next if $pos > $end;
		my $name=$a->qname;
		my $barcode=(split(/#/,$name))[-1];
		next if $barcode eq "0_0_0";
		$ref->{$barcode}=1;
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