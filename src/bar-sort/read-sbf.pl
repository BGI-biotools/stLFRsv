use strict;
use warnings;

open  IN,"$ARGV[0]" or die $!;
binmode(IN);
my $buf;
# seek(IN,405440,0);
while(read(IN,$buf,8)){
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
	
	my $bar1=$bar & 0xFFFFF;
	$bar =$bar >> 20;
	my $bar2=$bar & 0xFFFFF;
	$bar =$bar >> 20;
	my $bar3=$bar;
	
	# read(IN,$buf,8);
	# my $off=unpack("%64Q",$buf);
	
	print join(",",$bar3,$bar2,$bar1,$index,$name,$s,$e,$sup,"\n");
	# print join(",",$bar3,$bar2,$bar1,$off,"\n");
}
close(IN);
