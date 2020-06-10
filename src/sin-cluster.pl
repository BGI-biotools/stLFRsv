use strict;
use warnings;
use Cwd qw(abs_path getcwd);
use File::Basename;
use Statistics::R;

die "Usage: $0 <sbf file> <out-cluster file> <bin size>\n" if @ARGV != 3;
my ($seg,$out,$binsize)=@ARGV;


$out=abs_path($out);
my $dir=dirname($out);
my %hashL;
my %hashR;
open IN,"$seg" or die $!;
binmode(IN);
my $buf;
while(read(IN,$buf,8)){
	#read bin data
	my $bar=unpack("%64Q",$buf);

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
	######################
	
	my $posL=(int($s/$binsize)+1)*$binsize;
	my $posR=(int($e/$binsize))*$binsize;
	my $id=$bar."-".$index.":".$sup;
	$hashL{$name}{$posL}{"count"}++;
	push @{$hashL{$name}{$posL}{"id"}},$id;
	
	$hashR{$name}{$posR}{"count"}++;
	push @{$hashR{$name}{$posR}{"id"}},$id;
}

close IN;

my $li=0;
open RAW,">$out.raw" or die $!;
foreach my $i(
sort {
	my $na=$a;
	my $nb=$b;
	$na=~ s/^chr//;
	$nb=~ s/^chr//;
	$na=23 if $na eq "X";
	$na=24 if $na eq "Y";
	$na=25 if $na eq "M";
	$nb=23 if $nb eq "X";
	$nb=24 if $nb eq "Y";
	$nb=25 if $nb eq "M";
	if($na=~ /^\d+$/ and $nb=~ /^\d+$/){
		$na <=> $nb;
	}else{
		$na cmp $nb;
	}
}
keys %hashL){
	my $ni=0;
	foreach my $j(sort {$a <=> $b} keys %{$hashL{$i}}){
		my $line=join(",",@{$hashL{$i}{$j}{"id"}});
		print RAW "$li\t$ni\t$i\t$j\tL\t",$hashL{$i}{$j}{"count"},"\t$line\n";
		$ni++;
		$li++;
	}
}




foreach my $i(
sort {
	my $na=$a;
	my $nb=$b;
	$na=~ s/^chr//;
	$nb=~ s/^chr//;
	$na=23 if $na eq "X";
	$na=24 if $na eq "Y";
	$na=25 if $na eq "M";
	$nb=23 if $nb eq "X";
	$nb=24 if $nb eq "Y";
	$nb=25 if $nb eq "M";
	if($na=~ /^\d+$/ and $nb=~ /^\d+$/){
		$na <=> $nb;
	}else{
		$na cmp $nb;
	}
}
keys %hashR){
	my $ni=0;
	foreach my $j(sort {$a <=> $b} keys %{$hashR{$i}}){
		my $line=join(",",@{$hashR{$i}{$j}{"id"}});
		print RAW "$li\t$ni\t$i\t$j\tR\t",$hashR{$i}{$j}{"count"},"\t$line\n";
		$ni++;
		$li++;
	}
}
close RAW;

# %hashR=();
# %hashL=();

# if($low <= 0){
	# die "Wrong parameter for low depth filter!\n";
# }elsif($low < 1){
	# open IN,"$out.raw" or die $!;
	# open DAT,">$dir/sup.dat" or die $!;
	# while(<IN>){
		# chomp;
		# my @t=split;
		# print DAT "$t[5]\n";
	# }
	# close IN;
	# close DAT;
	# my $R = Statistics::R->new();
	# my $cmd=qq{
# library(MASS)
# a = read.table("$dir/sup.dat")
# x<-a[,1]
# subx<-x[which(x>(quantile(x,$ex)))]
# res=fitdistr(subx, "lognormal")
# th1<-qlnorm($low,res[[1]][1],res[[1]][2],lower.tail = TRUE)
# th2<-quantile(x,$low)
# };
	# $R->run($cmd);
	# my $low1=$R->get('th1');
	# my $low2=$R->get('th2');
	# $R->stop();
	# system("rm -rf $dir/sup.dat");
	# print STDERR "For single end at $low confidence:\n";
	# print STDERR "threshold1 is set to $low1\n";
	# print STDERR "threshold2 is set to $low2\n";
	# print STDERR "final threshold is set to ".int(($low1+$low2)/2 + 0.5)."\n";
	# if(abs($low1-$low2) >= 10){
		# print STDERR "Warning: too large diff between threshold1 and threshold2, the result may be unreliable!!!\n";
	# }
	# $low=int(($low1+$low2)/2 + 0.5);
# }

#get the distribu

# my @one_bat;
# my $lastc="";
# my $lasto="";
# my $lastp=-1;
# my $id=0;
# open OUT,">$out" or die $!;
# open IN,"$out.raw" or die $!;
# while(<IN>){
	# chomp;
	# my ($li,$ni,$chr,$pos,$ori,$count,$barline)=split;
	# if($count >= $low){
		# if(($chr ne $lastc) or ($ori ne $lasto)){
			# if(@one_bat > 0){
				# if(@one_bat < $nmerge +4){#too long outline bins
					# my ($chrf,$posf,$tcount,$tline,$ttline);
					# if($lasto eq "L"){
						# while(@one_bat > $nmerge){
							# pop @one_bat;
						# }
						# $chrf=$one_bat[0]->[0];
						# $posf=$one_bat[0]->[1];
						# $tcount=$one_bat[0]->[2];
						# $tline=$one_bat[0]->[3];
						# $ttline=$one_bat[0]->[4];
						# for(my $n=1;$n<=$#one_bat;$n++){
							# $tcount += $one_bat[$n]->[2];
							# $tline =$tline.",".$one_bat[$n]->[3];
							# $ttline =$ttline.",".$one_bat[$n]->[4];
						# }
					# }else{
						# while(@one_bat > $nmerge){
							# shift @one_bat;
						# }
						# $chrf=$one_bat[-1]->[0];
						# $posf=$one_bat[-1]->[1];
						# $tcount=$one_bat[-1]->[2];
						# $tline=$one_bat[-1]->[3];
						# $ttline=$one_bat[-1]->[4];
						# for(my $n=$#one_bat-1;$n>=0;$n--){
							# $tcount += $one_bat[$n]->[2];
							# $tline =$tline.",".$one_bat[$n]->[3];
							# $ttline =$ttline.",".$one_bat[$n]->[4];
						# }
					# }
					
					# print OUT "$id\t$chrf\t$posf\t$lasto\t$tcount\t$tline\t$ttline\n";
					# $id++;
				# }
			# }
			# @one_bat=();
			# push @one_bat,[($chr,$pos,$count,$ni,$li)];
			# $lastc=$chr;
			# $lasto=$ori;
			# $lastp=$pos;
		# }else{
			# if($pos-$lastp > $binsize){#allow one bin gap
				# if(@one_bat > 0){
					# if(@one_bat < $nmerge +4){#too long outline bins
						# my ($chrf,$posf,$tcount,$tline,$ttline);
						# if($lasto eq "L"){
							# while(@one_bat > $nmerge){
								# pop @one_bat;
							# }
							# $chrf=$one_bat[0]->[0];
							# $posf=$one_bat[0]->[1];
							# $tcount=$one_bat[0]->[2];
							# $tline=$one_bat[0]->[3];
							# $ttline=$one_bat[0]->[4];
							# for(my $n=1;$n<=$#one_bat;$n++){
								# $tcount += $one_bat[$n]->[2];
								# $tline =$tline.",".$one_bat[$n]->[3];
								# $ttline =$ttline.",".$one_bat[$n]->[4];
							# }
						# }else{
							# while(@one_bat > $nmerge){
								# shift @one_bat;
							# }
							# $chrf=$one_bat[-1]->[0];
							# $posf=$one_bat[-1]->[1];
							# $tcount=$one_bat[-1]->[2];
							# $tline=$one_bat[-1]->[3];
							# $ttline=$one_bat[-1]->[4];
							# for(my $n=$#one_bat-1;$n>=0;$n--){
								# $tcount += $one_bat[$n]->[2];
								# $tline =$tline.",".$one_bat[$n]->[3];
								# $ttline =$ttline.",".$one_bat[$n]->[4];
							# }
						# }
						
						# print OUT "$id\t$chrf\t$posf\t$lasto\t$tcount\t$tline\t$ttline\n";
						# $id++;
					# }
				# }
				# @one_bat=();
			# }
			# push @one_bat,[($chr,$pos,$count,$ni,$li)];
			# $lastp=$pos;	
		# }
	# }
# }

# if(@one_bat > 0){
	# if(@one_bat < $nmerge +4){#too long outline bins
		# my ($chrf,$posf,$tcount,$tline,$ttline);
		# if($lasto eq "L"){
			# while(@one_bat > $nmerge){
				# pop @one_bat;
			# }
			# $chrf=$one_bat[0]->[0];
			# $posf=$one_bat[0]->[1];
			# $tcount=$one_bat[0]->[2];
			# $tline=$one_bat[0]->[3];
			# $ttline=$one_bat[0]->[4];
			# for(my $n=1;$n<=$#one_bat;$n++){
				# $tcount += $one_bat[$n]->[2];
				# $tline =$tline.",".$one_bat[$n]->[3];
				# $ttline =$ttline.",".$one_bat[$n]->[4];
			# }
		# }else{
			# while(@one_bat > $nmerge){
				# shift @one_bat;
			# }
			# $chrf=$one_bat[-1]->[0];
			# $posf=$one_bat[-1]->[1];
			# $tcount=$one_bat[-1]->[2];
			# $tline=$one_bat[-1]->[3];
			# $ttline=$one_bat[-1]->[4];
			# for(my $n=$#one_bat-1;$n>=0;$n--){
				# $tcount += $one_bat[$n]->[2];
				# $tline =$tline.",".$one_bat[$n]->[3];
				# $ttline =$ttline.",".$one_bat[$n]->[4];
			# }
		# }
		
		# print OUT "$id\t$chrf\t$posf\t$lasto\t$tcount\t$tline\t$ttline\n";
		# $id++;
	# }
# }
# @one_bat=();

# close OUT;