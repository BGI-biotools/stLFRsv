use strict;
use warnings;
use Bio::DB::HTS;
use Statistics::Descriptive;
use Statistics::R;
use threads;
use threads::shared;
use Thread::Queue;

die "Usage: $0 <link all file> <split id file> <indexed bam file> <seg_len_freq> <phase info dir> <out file> <check size> <gap size> <merge_bin> <ncpu>\n" if @ARGV != 10;

my ($cluster,$split,$bamfile,$freqfile,$phase_dir,$outfile,$size,$gap,$Nmerge,$ncpu)=@ARGV;

my %reflen;
if(1){
	my $hfile = Bio::DB::HTSfile->open($bamfile);
	my $header = $hfile->header_read;
	my $n_targets = $header->n_targets;
	for(my $i=0;$i<$n_targets;$i++){
		my $name= $header->target_name->[$i];
		my $len= $header->target_len->[$i];
		$reflen{$name}=$len;
	}
}

my %seg_freq;
if(1){
	open IN,"$freqfile" or die $!;
	while(<IN>){
		chomp;
		my @t=split;
		$seg_freq{$t[0]}=$t[1];
	}
	close IN;
}



my %clusterinfo:shared;
open IN,"$cluster" or die $!;
while(<IN>){
	chomp;
	my @t=split;
	next unless $t[2] eq "PASS";
	$clusterinfo{"$t[0]-$t[1]"}=shared_clone([@t[3..8,11]]);
}
close IN;

my %svinfo;
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

my $TERM :shared = 0;
my $IDLE_QUEUE = Thread::Queue->new();
my $R_QUEUE = Thread::Queue->new();
$R_QUEUE->enqueue($1);
# $SIG{'INT'} = $SIG{'TERM'} =
# sub {
    # print(">>> Terminating <<<\n");
    # $TERM = 1;
    # $IDLE_QUEUE->insert(0, -1);
# };

if(1){
	my %work_queues;
	for (1..$ncpu) {
		my $work_q = Thread::Queue->new();
		my $thr = threads->create(\&worker, $work_q);
		$work_queues{$thr->tid()} = $work_q;
	}

	foreach my $svid(sort {$a <=> $b} keys %svinfo){
		my $tid;
		if(1){
			lock($IDLE_QUEUE);
			$tid= $IDLE_QUEUE->dequeue();
		}
		last if ($tid < 0);
		$work_queues{$tid}->enqueue($svid);
	}
	
	while(1){
		my $idle;
		if(1){
			lock($IDLE_QUEUE);
			$idle=$IDLE_QUEUE->pending();
		}
		if(defined $idle){
			last if $idle == $ncpu;
		}
		sleep(1);
	}
	
	$work_queues{$_}->enqueue(-1) foreach keys(%work_queues);
	
	$_->join() foreach threads->list();	
}

open OUT,">$outfile" or die $!;
foreach my $svid(sort {
	my @ta=split(/\t/,$svresult{$a});
	my @tb=split(/\t/,$svresult{$b});
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
}keys %svresult){
	print OUT $svresult{$svid}."\n";
}
close OUT;

sub worker{
    my ($work_q) = @_;
    my $tid = threads->tid();
	my $R=Statistics::R->new();
	my $hfile = Bio::DB::HTSfile->open($bamfile);
	my $index = Bio::DB::HTSfile->index_load($hfile);
	my $header = $hfile->header_read;
	my $stat = Statistics::Descriptive::Full->new();
    do {
		if(1){
			lock($IDLE_QUEUE);
			$IDLE_QUEUE->enqueue($tid);
		}
		my $svid = $work_q->dequeue();
		return if ($svid < 0);
		&processOne($svid,$R,$stat,$hfile,$index,$header);

    } while (! $TERM);

	return;
}

sub processOne{
	my ($svid,$R,$stat,$hfile,$index,$header)=@_;
	my $sv=$svinfo{$svid};
	
	my @result;
	&checkfour($sv,\@result,$R,$stat,$hfile,$index,$header);
	
	
	my $final;
	
	my ($score1,$score2)=(0,0);
	my ($cp,$cf)=(0,0);
	
	#mapQ
	if($result[2]=~ /PASS/ ){
		$score1 +=1;
		$cp++;
	}
	
	if($result[2]=~ /SUSPECT/ ){
		$score1 +=1;
		$cp++;
		$score2 +=(-1);
		$cf++;
	}
	
	if($result[2]=~ /FAILED/ ){
		$score2 +=(-1);
		$cf++;
	}
	#phase
	if($result[1]=~ /PASS/ ){
		$score1 +=1;
		$cp++;
	}
	
	if($result[1]=~ /NULL/ and $phase_dir ne "NULL"){
		$score2 +=(-1);
		$cf++;
	}
	
	if($result[1]=~ /FAILED/ ){
		$score2 +=(-2);
		$cf++;
	}
	
	my ($sym,$l1,$l2)=split(/,/,$result[4]);
	#sym
	if($sym=~ /SYM/){
		$score1 +=1;
		$cp++;
	}
	
	if($sym=~ /ASYM/){
		$score2 +=(-1);
		$cf++;
	}
	
	if($l1=~ /LONGER/){
		$score1 +=1;
		$cp++;
	}
	
	if($l1=~ /LOCAL/){
		$score2 +=(-1);
		$cf++;
	}
	
	if($l2=~ /LONGER/){
		$score1 +=1;
		$cp++;
	}
	
	if($l2=~ /LOCAL/){
		$score2 +=(-1);
		$cf++;
	}
	
	
	
	#PE
	if($result[3]=~ /PASS/ ){
		$score1 +=3;
		$cp++;
	}
	
	if($result[3]=~ /FAILED/ ){
		$score2 +=(-2);
		$cf++;
	}
	
	if($result[3]=~ /NULL/ ){
		$score2 +=(-1);
		$cf++;
	}
	
	
	#heat
	if($result[0]=~ /PASS/ ){
		$score1 +=4;
		$cp++;
	}
	
	if($result[0]=~ /SUSPECT/ ){
		$score1 +=3;
		$cp++;
	}
	
	if($result[0]=~ /LOW/ ){
		$score2 +=(-2);
		$cf++;
	}
	
	my $score=$score1*$cp+$score2*$cf;
	
	if($phase_dir ne "NULL"){
		if($score >=35 ){
			$final="PASS";
		}else{
			$final="*";
		}	
	}else{
		if($score >=20 ){
			$final="PASS";
		}else{
			$final="*";
		}
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
	my $line=join("\t",$svid,@{$sv}[0..6],$simple_type,$score,$final,@result);
	
	if(1){
		lock(%svresult);
		$svresult{$svid}=$line;
	}

	return;
}

sub checkfour{
	my ($ref,$out,$R,$stat,$hfile,$index,$header)=@_;
	
	##1 heatmap
	my (@dep1,@dep2);
	my ($start,$end);
	my ($type1,$type2)=split(//,$ref->[6]);
	
	if($type1 eq "L"){
		$start=$ref->[3];
		$start=0 if $start <0;
		$end=$ref->[3]+10*$size;
		$end=$reflen{$ref->[2]} if $end > $reflen{$ref->[2]};
	}else{
		$start=$ref->[3]-10*$size;
		$start=0 if $start <0;
		$end=$ref->[3];
		$end=$reflen{$ref->[2]} if $end > $reflen{$ref->[2]};
	}

	my @bararray;
	&get_bar_from_bam($ref->[2],$start,$end,\@bararray,$hfile,$index,$header);
	foreach my $r(@bararray){
		my $count= keys %{$r};
		push @dep1,$count;
	}
	
	@bararray=();
	
	if($type2 eq "L"){
		$start=$ref->[5];
		$start=0 if $start <0;
		$end=$ref->[5]+10*$size;
		$end=$reflen{$ref->[4]} if $end > $reflen{$ref->[4]};
	}else{
		$start=$ref->[5]-10*$size;
		$start=0 if $start <0;
		$end=$ref->[5];
		$end=$reflen{$ref->[4]} if $end > $reflen{$ref->[4]};
	}

	&get_bar_from_bam($ref->[4],$start,$end,\@bararray,$hfile,$index,$header);
	foreach my $r(@bararray){
		my $count= keys %{$r};
		push @dep2,$count;
	}
	
	my $judge1="";
	unless (@dep1 > 5 and @dep2 > 5){
		$judge1 = "NULL";
		goto PASS;
	}
	
	$stat->clear();
	$stat->add_data(@dep1);
	my $avg_bar_dep1=$stat->mean();
	
	$stat->clear();
	$stat->add_data(@dep2);
	my $avg_bar_dep2=$stat->mean();
	
	##freq
	my($max_len1,$max_len2)=(0,0);
	my (@exp1,@exp2);
	foreach my $i(sort {$a <=> $b} keys %seg_freq){
		my $exp=$seg_freq{$i}*$avg_bar_dep1/2;
		if($exp >= 4){
			$max_len1=$i;
			push @exp2,$exp;
		}else{
			last;
		}
	}
	
	foreach my $i(sort {$a <=> $b} keys %seg_freq){
		my $exp=$seg_freq{$i}*$avg_bar_dep2/2;
		if($exp >= 4){
			$max_len2=$i;
			push @exp1,$exp;
		}else{
			last;
		}
	}
	
	my $max_len= $max_len1 < $max_len2 ? $max_len1 :$max_len2;
	
	unless($max_len >= 5*$size){
		$judge1 = "NULL";
		goto PASS;
	}
	
	##
	my $max_bin=int($max_len/$size);
	my @matrix;
	my (@bar1,@bar2);
	my $f_bin=&get_matrix($ref->[2],$ref->[3],$ref->[4],$ref->[5],$max_bin,$size,\@bar1,\@bar2,$hfile,$index,$header);
	
	unless ($f_bin >=5){
		$judge1 = "NULL";
		goto PASS;
	}
	###
	my ($r_bin,$mode);
	if($ref->[2] eq $ref->[4]){	
		$r_bin= $f_bin < int(($ref->[5]-$ref->[3])/$size/4) ? $f_bin : int(($ref->[5]-$ref->[3])/$size/4);
		if($r_bin <5){
			$r_bin=5;
			$mode="C";
		}else{
			$mode="A";
		}
		$r_bin =30 if $r_bin > 30;
	}else{
		$mode="A";
		$r_bin=$f_bin < int($gap/$size+0.5) ? $f_bin : int($gap/$size+0.5);
		$r_bin =30 if $r_bin > 30;
		$r_bin =5 if $r_bin < 5;
	}
	
	#share
	my %barhash;
	my %tmphash;
	for(my $i=$f_bin- $r_bin;$i<= $f_bin+$r_bin;$i++){
		foreach my $bar(keys %{$bar1[$i]}){
			push @{$barhash{$bar}},$i;
		}
		
		foreach my $bar(keys %{$bar2[$i]}){
			push @{$barhash{$bar}},$i+2*$f_bin+1;
		}
	
	}
	
	for(my $i=$f_bin- $r_bin;$i<= $f_bin+$r_bin;$i++){
		%tmphash=();
		foreach my $bar(keys %{$bar1[$i]}){
			foreach my $j (@{$barhash{$bar}}){
				next unless $i < $j;
				$tmphash{$j}++;
			}
		}
		
		for(my $j=$i+1;$j<=$f_bin+$r_bin;$j++){
			if(exists $tmphash{$j}){
				$matrix[$i][$j]=$tmphash{$j};
			}else{
				$matrix[$i][$j]=0;
			}
		}
		
		for(my $j=$f_bin- $r_bin;$j<= $f_bin+$r_bin;$j++){
			my $j2=$j+2*$f_bin+1;
			if(exists $tmphash{$j2}){
				$matrix[$i][$j2]=$tmphash{$j2};
			}else{
				$matrix[$i][$j2]=0;
			}
		}
		
		%tmphash=();
		foreach my $bar(keys %{$bar2[$i]}){
			foreach my $j (@{$barhash{$bar}}){
				next unless ($i+2*$f_bin+1) < $j;
				$tmphash{$j}++;
			}
		}
		
		for(my $j=$i+1;$j<=$f_bin+$r_bin;$j++){
			my $j2=$j+2*$f_bin+1;
			if(exists $tmphash{$j2}){
				$matrix[$i+2*$f_bin+1][$j2]=$tmphash{$j2};
			}else{
				$matrix[$i+2*$f_bin+1][$j2]=0;
			}
		}
	}
	#####
	
	my (@r1,@r2,@r3,@r4);
	for(my $i=0;$i<$r_bin;$i++){
		my $off1=$f_bin-1;
		for(my $j=0;$j<$r_bin;$j++){
			my $off2=3*$f_bin+2;
			push @r1,$matrix[$off1-$i][$off2+$j];
		}
	}
	
	for(my $i=0;$i<$r_bin;$i++){
		my $off1=$f_bin+1;
		for(my $j=0;$j<$r_bin;$j++){
			my $off2=3*$f_bin+2;
			push @r2,$matrix[$off1+$i][$off2+$j];
		}
	}
	
	for(my $i=0;$i<$r_bin;$i++){
		my $off1=$f_bin+1;
		for(my $j=0;$j<$r_bin;$j++){
			my $off2=3*$f_bin;
			push @r3,$matrix[$off1+$i][$off2-$j];
		}
	}
	
	for(my $i=0;$i<$r_bin;$i++){
		my $off1=$f_bin-1;
		for(my $j=0;$j<$r_bin;$j++){
			my $off2=3*$f_bin;
			push @r4,$matrix[$off1-$i][$off2-$j];
		}
	}
	
	my (@r13,@r24);
	push @r13,@r1;
	push @r13,@r3;
	push @r24,@r2;
	push @r24,@r4;
	
	my (@ex1,@ex2);
	#1
	if($type1 eq "L"){
		%barhash=();
		for(my $i=0;$i<$f_bin;$i++){
			foreach my $bar(keys %{$bar1[$f_bin+1+$i]}){
				push @{$barhash{$bar}},$i;
			}
		}
		
		%tmphash=();
		foreach my $bar(keys %{$bar2[$f_bin]}){
			foreach my $i (@{$barhash{$bar}}){
				$tmphash{$i}++;
			}
		}
		
		for(my $i=0;$i<$f_bin;$i++){
			if(exists $tmphash{$i}){
				push @ex1,$tmphash{$i};
			}else{
				push @ex1,0;
			}
		}
	}else{
		%barhash=();
		for(my $i=$f_bin-1;$i>=0;$i--){
			foreach my $bar(keys %{$bar1[$i]}){
				push @{$barhash{$bar}},$i;
			}
		}
		
		%tmphash=();
		foreach my $bar(keys %{$bar2[$f_bin]}){
			foreach my $i (@{$barhash{$bar}}){
				$tmphash{$i}++;
			}
		}
		
		for(my $i=$f_bin-1;$i>=0;$i--){
			if(exists $tmphash{$i}){
				push @ex1,$tmphash{$i};
			}else{
				push @ex1,0;
			}
		}
	}
	
	#2
	if($type2 eq "L"){
		%barhash=();
		for(my $i=0;$i<$f_bin;$i++){
			foreach my $bar(keys %{$bar2[$f_bin+1+$i]}){
				push @{$barhash{$bar}},$i;
			}
		}
		
		%tmphash=();
		foreach my $bar(keys %{$bar1[$f_bin]}){
			foreach my $i (@{$barhash{$bar}}){
				$tmphash{$i}++;
			}
		}
		
		for(my $i=0;$i<$f_bin;$i++){
			if(exists $tmphash{$i}){
				push @ex2,$tmphash{$i};
			}else{
				push @ex2,0;
			}
		}
	}else{
		%barhash=();
		for(my $i=$f_bin-1;$i>=0;$i--){
			foreach my $bar(keys %{$bar2[$i]}){
				push @{$barhash{$bar}},$i;
			}
		}
		
		%tmphash=();
		foreach my $bar(keys %{$bar1[$f_bin]}){
			foreach my $i (@{$barhash{$bar}}){
				$tmphash{$i}++;
			}
		}
		
		for(my $i=$f_bin-1;$i>=0;$i--){
			if(exists $tmphash{$i}){
				push @ex2,$tmphash{$i};
			}else{
				push @ex2,0;
			}
		}
	}
	
	my $cmd;
	$cmd=qq{rm(list = ls())
r1<-c()
r2<-c()
r3<-c()
r4<-c()
cr1<-c()
cr2<-c()
ex1<-c()
ex2<-c()
};
	
	for(my $i=0;$i<=$#r1;$i+=100){
		my $start=$i;
		my $end=$i+99;
		$end = $#r1 if $end > $#r1;
		my $line=join(",",@r1[$start..$end]);
		$cmd.="r1<-c(r1,$line)\n";
	}
	
	for(my $i=0;$i<=$#r2;$i+=100){
		my $start=$i;
		my $end=$i+99;
		$end = $#r2 if $end > $#r2;
		my $line=join(",",@r2[$start..$end]);
		$cmd.="r2<-c(r2,$line)\n";
	}
	
	for(my $i=0;$i<=$#r3;$i+=100){
		my $start=$i;
		my $end=$i+99;
		$end = $#r3 if $end > $#r3;
		my $line=join(",",@r3[$start..$end]);
		$cmd.="r3<-c(r3,$line)\n";
	}
	
	for(my $i=0;$i<=$#r4;$i+=100){
		my $start=$i;
		my $end=$i+99;
		$end = $#r4 if $end > $#r4;
		my $line=join(",",@r4[$start..$end]);
		$cmd.="r4<-c(r4,$line)\n";
	}
	
	for(my $i=0;$i<=$#r13;$i+=100){
		my $start=$i;
		my $end=$i+99;
		$end = $#r13 if $end > $#r13;
		my $line=join(",",@r13[$start..$end]);
		$cmd.="cr1<-c(cr1,$line)\n";
	}
	
	for(my $i=0;$i<=$#r24;$i+=100){
		my $start=$i;
		my $end=$i+99;
		$end = $#r24 if $end > $#r24;
		my $line=join(",",@r24[$start..$end]);
		$cmd.="cr2<-c(cr2,$line)\n";
	}
	
	for(my $i=0;$i<=$#ex1;$i+=100){
		my $start=$i;
		my $end=$i+99;
		$end = $#ex1 if $end > $#ex1;
		my $line=join(",",@ex1[$start..$end]);
		$cmd.="ex1<-c(ex1,$line)\n";
	}
	
	for(my $i=0;$i<=$#ex2;$i+=100){
		my $start=$i;
		my $end=$i+99;
		$end = $#ex2 if $end > $#ex2;
		my $line=join(",",@ex2[$start..$end]);
		$cmd.="ex2<-c(ex2,$line)\n";
	}

	$cmd.=qq{h12<-wilcox.test(r1,r2,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
h13<-wilcox.test(r1,r3,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
h14<-wilcox.test(r1,r4,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
h23<-wilcox.test(r2,r3,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
h24<-wilcox.test(r2,r4,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
h34<-wilcox.test(r3,r4,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
h_two<-wilcox.test(cr1,cr2,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
h_ex<-wilcox.test(ex1,ex2,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
};	
	my ($p12,$p13,$p14,$p23,$p24,$p34,$p_two,$p_ex);
	if(1){
		lock($R_QUEUE);
		$R_QUEUE->dequeue();
		$R->run($cmd);
		$p12=$R->get('h12$p.value');
		$p13=$R->get('h13$p.value');
		$p14=$R->get('h14$p.value');
		$p23=$R->get('h23$p.value');
		$p24=$R->get('h24$p.value');
		$p34=$R->get('h34$p.value');
		$p_two=$R->get('h_two$p.value');
		$p_ex=$R->get('h_ex$p.value');
		$R_QUEUE->enqueue(1);
	}
	
	
	my $type=$ref->[6];
	my $score=0;
	if($type eq "RL"){
		my ($s1,$s2)=(0,0);
		my ($ck,$cf);
		if($mode eq "C"){
			#T1
			$ck=0;
			$cf=0;
			
			if($p12 <= 1e-5){
				$s1+=2.5*($ck+1);
				$ck++;
			}elsif($p12 <=1e-2){
				$s1+=1.5*($ck+1);
				$ck++;
			}elsif($p12 > 0.5){
				$s1-=1.5*($cf+1);
				$cf++;
			}
			
			
			if($p14 <= 1e-5){
				$s1+=2.5*($ck+1);
				$ck++;
			}elsif($p14 <=1e-2){
				$s1+=1.5*($ck+1);
				$ck++;
			}elsif($p14 > 0.5){
				$s1-=1.5*($cf+1);
				$cf++;
			}
			
			if($p13 <= 1e-5){
				$s1+=3*($ck+1);
				$ck++;
			}elsif($p13 <=1e-2){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p13 > 0.5){
				$s1-=1*($cf+1);
				$cf++;
			}
			
			$s1=int($s1+0.5);
		}else{
			#T1
			$ck=0;
			$cf=0;
			if($p12 <= 1e-5){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p12 <=1e-2){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p12 > 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p13 <= 1e-5){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p13 <=1e-2){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p13 > 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p14 <= 1e-5){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p14 <=1e-2){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p14 > 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			#T2
			$ck=0;
			$cf=0;
			
			if($p23 >=(1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p23 < 0.5){
				$s2-=1*($cf+1);
				$cf++;
			}
			
			if($p34 <= 1e-5){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p34 > 0.5){
				$s2-=1*($cf+1);
				$cf++;
			}
			
			if($p13 > 1e-2 and $p13 < (1-1e-2)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p13 > 1e-5 and $p13 < (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p13 <= 1e-5 or $p13 >= (1-1e-5)){
				$s2-=2*($cf+1);
				$cf++;
			}
			
			if($p_two <= 1e-5){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p_two <= 1e-2){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p_two > 0.5){
				$s2-=2*($cf+1);
				$cf++;
			}
			
		}
	
		$score= $s1 > $s2 ? $s1 : $s2;
	}elsif($type eq "LR"){
		my ($s1,$s2)=(0,0);
		my $ck;
		my $cf;
		if($mode eq "C"){
			#T1
			$ck=0;
			$cf=0;
			
			if($p13 >= (1-1e-5)){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p13 >=(1-1e-2)){
				$s1+=0.5*($ck+1);
				$ck++;
			}elsif($p13 < 0.5){
				$s1-=3*($cf+1);
				$cf++;
			}
			
			if($p23 >= (1-1e-5)){
				$s1+=1.5*($ck+1);
				$ck++;
			}elsif($p23 >=(1-1e-2)){
				$s1+=0.5*($ck+1);
				$ck++;
			}elsif($p23 < 0.5){
				$s1-=2.5*($cf+1);
				$cf++;
			}
			
			if($p34 <= 1e-5){
				$s1+=1.5*($ck+1);
				$ck++;
			}elsif($p34 <= 1e-2){
				$s1+=0.5*($ck+1);
				$ck++;
			}elsif($p34 > 0.5){
				$s1-=2.5*($cf+1);
				$cf++;
			}
			
			$s1=int($s1+0.5);
		}else{
			#T1
			$ck=0;
			$cf=0;
			if($p13 >= (1-1e-5)){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p13 >=(1-1e-2)){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p13 < 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p23 >= (1-1e-5)){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p23 >= (1-1e-2)){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p23 < 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p34 <= 1e-5){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p34 <= 1e-2){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p34 > 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			#T2
			$ck=0;
			$cf=0;
			
			if($p12 <= 1e-5){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p12 > 0.5){
				$s2-=1*($cf+1);
				$cf++;
			}
			
			if($p14 <= 1e-5){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p14 > 0.5){
				$s2-=1*($cf+1);
				$cf++;
			}
			
			if($p13 > 1e-2 and $p13 < (1-1e-2)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p13 > 1e-5 and $p13 < (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p13 <= 1e-5 or $p13 >= (1-1e-5)){
				$s2-=2*($cf+1);
				$cf++;
			}
			
			if($p_two <= 1e-5){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p_two <= 1e-2){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p_two > 0.5){
				$s2-=2*($cf+1);
				$cf++;
			}
			
		}
		
		$score= $s1 > $s2 ? $s1 : $s2;
	}elsif($type eq "LL"){
		my ($s1,$s2)=(0,0);
		my $ck;
		my $cf;
		if($mode eq "C"){
			#T1
			$ck=0;
			$cf=0;
			if($p12 >= (1-1e-5)){
				$s1+=1.5*($ck+1);
				$ck++;
			}elsif($p12 >=(1-1e-2)){
				$s1+=0.5*($ck+1);
				$ck++;
			}elsif($p12 < 0.5){
				$s1-=2.5*($cf+1);
				$cf++;
			}
			
			if($p24 <= 1e-5){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p24 <= 1e-2){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p24 > 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p23 <= 1e-5){
				$s1+=2.5*($ck+1);
				$ck++;
			}elsif($p23 <= 1e-2){
				$s1+=1.5*($ck+1);
				$ck++;
			}elsif($p23 > 0.5){
				$s1-=1.5*($cf+1);
				$cf++;
			}
			
			#T2
			$ck=0;
			$cf=0;
			
			if($p14 >= (1-1e-5)){
				$s2+=0.5*($ck+1);
				$ck++;
			}elsif($p14 < 0.5){
				$s2-=1.5*($cf+1);
				$cf++;
			}
			
			if($p34 >= (1-1e-5)){
				$s2+=1.5*($ck+1);
				$ck++;
			}elsif($p34 < 0.5){
				$s2-=0.5*($cf+1);
				$cf++;
			}
			
			if($p24 > 1e-2 and $p24 < (1-1e-2)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p24 > 1e-5 and $p24 < (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p24 <= 1e-5 or $p24 >= (1-1e-5)){
				$s2-=2*($cf+1);
				$cf++;
			}
			
			
			if($p_two >= (1-1e-5)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p_two >= (1-1e-2)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p_two < 0.5){
				$s2-=2*($cf+1);
				$cf++;
			}
			
			
			$s1=int($s1+0.5);
			$s2=int($s2+0.5);
		}else{
			#T1
			$ck=0;
			$cf=0;
			if($p12 >= (1-1e-5)){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p12 >=(1-1e-2)){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p12 < 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p23 <= 1e-5){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p23 <= 1e-2){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p23 > 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p24 <= 1e-5){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p24 <=1e-2){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p24 > 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			#T2
			$ck=0;
			$cf=0;
			
			if($p14 >= (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p14 < 0.5){
				$s2-=1*($cf+1);
				$cf++;
			}
			
			if($p34 >= (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p34 < 0.5){
				$s2-=1*($cf+1);
				$cf++;
			}
			
			if($p24 > 1e-2 and $p24 < (1-1e-2)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p24 > 1e-5 and $p24 < (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p24 <= 1e-5 or $p24 >=(1-1e-5)){
				$s2-=2*($cf+1);
				$cf++;
			}
			
			if($p_two >= (1-1e-5)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p_two >= (1-1e-2)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p_two < 0.5){
				$s2-=2*($cf+1);
				$cf++;
			}
		}
	
		$score= $s1 > $s2 ? $s1 : $s2;
	}elsif($type eq "RR"){
		my ($s1,$s2)=(0,0);
		my $ck;
		my $cf;
		if($mode eq "C"){
			#T1
			$ck=0;
			$cf=0;
			if($p14 >= (1-1e-5)){
				$s1+=1.5*($ck+1);
				$ck++;
			}elsif($p14 >= (1-1e-2)){
				$s1+=0.5*($ck+1);
				$ck++;
			}elsif($p14 < 0.5){
				$s1-=2.5*($cf+1);
				$cf++;
			}
			
			if($p24 >= (1-1e-5)){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p24 >= (1-1e-2)){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p24 < 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p34 >= (1-1e-5)){
				$s1+=2.5*($ck+1);
				$ck++;
			}elsif($p34 >=(1-1e-2)){
				$s1+=1.5*($ck+1);
				$ck++;
			}elsif($p34 < 0.5){
				$s1-=1.5*($cf+1);
				$cf++;
			}
			
			#T2
			$ck=0;
			$cf=0;
			
			if($p12 >= (1-1e-5)){
				$s2+=0.5*($ck+1);
				$ck++;
			}elsif($p12 < 0.5){
				$s2-=1.5*($cf+1);
				$cf++;
			}
			
			if($p23 <= 1e-5){
				$s2+=1.5*($ck+1);
				$ck++;
			}elsif($p23 > 0.5){
				$s2-=0.5*($cf+1);
				$cf++;
			}
			
			if($p24 > 1e-2 and $p24 < (1-1e-2)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p24 > 1e-5 and $p24 < (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p24 <= 1e-5 or $p24 >= (1-1e-5)){
				$s2-=2*($cf+1);
				$cf++;
			}
			
			if($p_two >= (1-1e-5)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p_two >= (1-1e-2)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p_two < 0.5){
				$s2-=2*($cf+1);
				$cf++;
			}
			
			
			$s1=int($s1+0.5);
			$s2=int($s2+0.5);
		}else{
			#T1
			$ck=0;
			$cf=0;
			if($p14 >= (1-1e-5)){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p14 >=(1-1e-2)){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p14 < 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p24 >= (1-1e-5)){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p24 >= (1-1e-2)){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p24 < 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			if($p34 >= (1-1e-5)){
				$s1+=2*($ck+1);
				$ck++;
			}elsif($p34 >= (1-1e-2)){
				$s1+=1*($ck+1);
				$ck++;
			}elsif($p34 < 0.5){
				$s1-=2*($cf+1);
				$cf++;
			}
			
			#T2
			$ck=0;
			$cf=0;
			if($p12 >= (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p12 < 0.5){
				$s2-=1*($cf+1);
				$cf++;
			}
			
			if($p23 <= 1e-5){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p23 > 0.5){
				$s2-=1*($cf+1);
				$cf++;
			}
			
			if($p24 > 1e-2 and $p24 < (1-1e-2)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p24 > 1e-5 and $p24 < (1-1e-5)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p24 <= 1e-5 or $p24 >=(1-1e-5)){
				$s2-=2*($cf+1);
				$cf++;
			}
			
			if($p_two >= (1-1e-5)){
				$s2+=2*($ck+1);
				$ck++;
			}elsif($p_two >= (1-1e-2)){
				$s2+=1*($ck+1);
				$ck++;
			}elsif($p_two < 0.5){
				$s2-=2*($cf+1);
				$cf++;
			}
			
		}
		
		$score= $s1 > $s2 ? $s1 : $s2;
	}
	
	
	PASS:
	my $ex_line;
	if($judge1 eq "NULL"){
		$judge1="$judge1:NA";
		$ex_line="NULL,NA,NA";
	}else{
		if($score >10){
			$judge1="PASS:$score";
		}elsif($score >= 6){
			$judge1="SUSPECT:$score";
		}else{
			$judge1="LOW:$score";
		}
		
		my (@dex1,@dex2);
		for(my $i=0;$i<=$#ex1;$i++){
			my $delta=$ex1[$i]-($exp1[$i]);
			push @dex1,$delta;
		}
		
		for(my $i=0;$i<=$#ex2;$i++){
			my $delta=$ex2[$i]-($exp2[$i]);
			push @dex2,$delta;
		}
		
		my $len1=@dex1;
		my $len2=@dex2;
		
		$cmd=qq{
dex1<-c()
dex2<-c()
cex1<-rep(0,$len1)
cex2<-rep(0,$len2)
};
		for(my $i=0;$i<=$#dex1;$i+=100){
			my $start=$i;
			my $end=$i+99;
			$end = $#dex1 if $end > $#dex1;
			my $line=join(",",@dex1[$start..$end]);
			$cmd.="dex1<-c(dex1,$line)\n";
		}
		
		for(my $i=0;$i<=$#dex2;$i+=100){
			my $start=$i;
			my $end=$i+99;
			$end = $#dex2 if $end > $#dex2;
			my $line=join(",",@dex2[$start..$end]);
			$cmd.="dex2<-c(dex2,$line)\n";
		}
		
		$cmd.=qq{
h_dx1<-wilcox.test(dex1,cex1,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
h_dx2<-wilcox.test(dex2,cex2,paired=FALSE, conf.level = 0.95,alternative='g',exact=T,correct=F)
};
		my ($p_dx1,$p_dx2);
		if(1){
			lock($R_QUEUE);
			$R_QUEUE->dequeue();
			$R->run($cmd);
			$p_dx1=$R->get('h_dx1$p.value');
			$p_dx2=$R->get('h_dx2$p.value');
			$R_QUEUE->enqueue(1);
		}
		
		
		if($p_ex <= 1e-3 or $p_ex >=(1-1e-3)){
			$ex_line="ASYM"
		}else{
			$ex_line="SYM"
		}
		
		if($p_dx1 <= 0.05){
			$ex_line.=",LONGER:".($size/2+$f_bin*$size);
		}else{
			my $l_bin=&search_extend(\@dex1);
			if($type1 eq "L"){
				$ex_line.=",LOCAL1:$ref->[2]-".$ref->[3]."_".($ref->[3]+($size/2+$l_bin*$size));
			}else{
				$ex_line.=",LOCAL1:$ref->[2]-".($ref->[3]-($size/2+$l_bin*$size))."_".$ref->[3];
			}
		}
		
		if($p_dx2 <= 0.05){
			$ex_line.=",LONGER:".($size/2+$f_bin*$size);
		}else{
			my $l_bin=&search_extend(\@dex2);
			if($type2 eq "L"){
				$ex_line.=",LOCAL2:$ref->[4]-".$ref->[5]."_".($ref->[5]+($size/2+$l_bin*$size));
			}else{
				$ex_line.=",LOCAL2:$ref->[4]-".($ref->[5]-($size/2+$l_bin*$size))."_".$ref->[5];
			}
		}
	}

	
	$out->[0]=$judge1;
	
	#barcode phase HP
	my $judge2;
	my $chr1=$ref->[2];
	my $pos1=$ref->[3];
	my $check=1;
	
	my $hp1;
	if(exists $phaseregion{$chr1}){
		foreach my $i(keys %{$phaseregion{$chr1}}){
			if($pos1 >= $phaseregion{$chr1}{$i}->[0]  and $pos1 <= $phaseregion{$chr1}{$i}->[1]){
				$check=0;
				$hp1="${chr1}_${i}_".$phaseregion{$chr1}{$i}->[0]."-".$phaseregion{$chr1}{$i}->[1];
				last;
			}
		}
	}
	
	if($check){
		$judge2="NULL:REGION:NA:NA:NA:NA";
		goto PHASE;
	}
	
	$check=1;
	my $chr2=$ref->[4];
	my $pos2=$ref->[5];
	my $hp2;
	if(exists $phaseregion{$chr2}){
		foreach my $i(keys %{$phaseregion{$chr2}}){
			if($pos2 >= $phaseregion{$chr2}{$i}->[0]  and $pos2 <= $phaseregion{$chr2}{$i}->[1]){
				$check=0;
				$hp2="${chr2}_${i}_".$phaseregion{$chr2}{$i}->[0]."-".$phaseregion{$chr2}{$i}->[1];
				last;
			}	
		}
	}
	
	if($check){
		$judge2="NULL:REGION:NA:NA:NA:NA";
		goto PHASE;
	}
	
	unless($f_bin=~ /^\d+$/){
		$judge2="NULL:REGION:NA:NA:NA:NA";
		goto PHASE;
	}
	####
	
	my @share;
	if(exists $clusterinfo{$ref->[0]."-".$ref->[1]}){
		foreach my $i(split(/,/,$clusterinfo{$ref->[0]."-".$ref->[1]}->[6])){
			my $bar1=$i & 0xFFFFF;
			$i =$i >> 20;
			my $bar2=$i & 0xFFFFF;
			$i =$i >> 20;
			my $bar3=$i & 0xFFFFF;
			
			my $line=join("_",$bar3,$bar2,$bar1);
			
			push @share,$line;
		}
	}elsif(exists $clusterinfo{$ref->[1]."-".$ref->[0]}){
		foreach my $i(split(/,/,$clusterinfo{$ref->[1]."-".$ref->[0]}->[6])){
			my $bar1=$i & 0xFFFFF;
			$i =$i >> 20;
			my $bar2=$i & 0xFFFFF;
			$i =$i >> 20;
			my $bar3=$i & 0xFFFFF;
			
			my $line=join("_",$bar3,$bar2,$bar1);
			
			push @share,$line;
		}
	}else{
		die "Wrong linked all file!!!\n";
	}

	unless(@share >5){
		$judge2="NULL:COUNT:NA:NA:NA:NA";
		goto PHASE;
	}
	###
	
	my @stat1=(0,0,0,0);
	my @stat2=(0,0,0,0);
	foreach my $bar(@share){
		$stat1[0]+=1;
		$stat2[0]+=1;
		if($chr1 eq $chr2){
			if(exists $phaseinfo{$chr1}){
				if(exists $phaseinfo{$chr1}{$bar}){
					if($phaseinfo{$chr1}{$bar} eq "PAT"){
						$stat1[1]+=1;
						$stat2[1]+=1;
					}else{
						$stat1[2]+=1;
						$stat2[2]+=1;
					}
				}else{
					$stat1[3]+=1;
					$stat2[3]+=1;
				}
			}else{
				$stat1[3]+=1;
				$stat2[3]+=1;
			}
		}else{
			if(exists $phaseinfo{$chr1}){
				if(exists $phaseinfo{$chr1}{$bar}){
					if($phaseinfo{$chr1}{$bar} eq "PAT"){
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
			
			if(exists $phaseinfo{$chr2}){
				if(exists $phaseinfo{$chr2}{$bar}){
					if($phaseinfo{$chr2}{$bar} eq "PAT"){
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
	}
	
	my $ratio1=$stat1[3]/$stat1[0];
	my $ratio2=$stat2[3]/$stat2[0];
	if($ratio1 > 0.75){
		$judge2="NULL:UNPHASED:$hp1:$hp2:NA:NA";
		goto PHASE;
	}
	
	if($ratio2 > 0.75){
		$judge2="NULL:UNPHASED:$hp1:$hp2:NA:NA";
		goto PHASE;
	}
	
	if(($stat1[1]+$stat1[2])< 50){
		my $tnum=$stat1[1]+$stat1[2];
		$stat1[1]=int($stat1[1]*(50/$tnum));
		$stat1[2]=int($stat1[2]*(50/$tnum));
	}
	
	if(($stat2[1]+$stat2[2])< 50){
		my $tnum=$stat2[1]+$stat2[2];
		$stat2[1]=int($stat2[1]*(50/$tnum));
		$stat2[2]=int($stat2[2]*(50/$tnum));
	}
	
	my $mean1=int(($stat1[1]+$stat1[2])/2);
	my $mean2=int(($stat2[1]+$stat2[2])/2);
	my $tnum1=$stat1[1]+$stat1[2];
	my $tnum2=$stat2[1]+$stat2[2];
	
	$cmd=qq{
a1<-fisher.test(matrix(c($stat1[1],$stat1[2],$tnum1,0),nrow=2),hybrid = FALSE,alternative = "t")
b1<-fisher.test(matrix(c($stat1[1],$stat1[2],0,$tnum1),nrow=2),hybrid = FALSE,alternative = "t")
c1<-fisher.test(matrix(c($stat1[1],$stat1[2],$mean1,$mean1),nrow=2),hybrid = FALSE,alternative = "t")
a2<-fisher.test(matrix(c($stat2[1],$stat2[2],$tnum2,0),nrow=2),hybrid = FALSE,alternative = "t")
b2<-fisher.test(matrix(c($stat2[1],$stat2[2],0,$tnum2),nrow=2),hybrid = FALSE,alternative = "t")
c2<-fisher.test(matrix(c($stat2[1],$stat2[2],$mean2,$mean2),nrow=2),hybrid = FALSE,alternative = "t")
};	
	my($pa1,$pb1,$pc1,$pa2,$pb2,$pc2);
	if(1){
		lock($R_QUEUE);
		$R_QUEUE->dequeue();
		$R->run($cmd);
		$pa1=$R->get('a1$p.value');
		$pb1=$R->get('b1$p.value');
		$pc1=$R->get('c1$p.value');
		$pa2=$R->get('a2$p.value');
		$pb2=$R->get('b2$p.value');
		$pc2=$R->get('c2$p.value');
		$R_QUEUE->enqueue(1);
	}
	
	
	my (@ty1,@ty2);
	unless($pa1 < 1e-5){
		push @ty1,"1|0";
	}
	
	unless($pb1 < 1e-5){
		push @ty1,"0|1";
	}
	
	unless($pc1 < 0.05){
		push @ty1,"1|1";
	}
	
	unless($pa2 < 1e-5){
		push @ty2,"1|0";
	}
	
	unless($pb2 < 1e-5){
		push @ty2,"0|1";
	}
	
	unless($pc2 < 0.05){
		push @ty2,"1|1";
	}
	
	if(@ty1 != 1){
		$ty1[0]="FAILED";
	}
	
	if(@ty2 != 1){
		$ty2[0]="FAILED";
	}
	
	if($chr1 eq $chr2){
		if($ty1[0] eq "FAILED"){
			$judge2="FAILED:MIX:$hp1:$hp2:$pa1,$pb1,$pc1:$pa2,$pb2,$pc2";
		}else{
			$judge2="PASS:($ty1[0])-($ty2[0]):$hp1:$hp2:$pa1,$pb1,$pc1:$pa2,$pb2,$pc2";
		}
	}else{
		if($ty1[0] eq $ty2[0]){
			if($ty1[0] eq "FAILED"){
				$judge2="FAILED:MIX:$hp1:$hp2:$pa1,$pb1,$pc1:$pa2,$pb2,$pc2";
			}else{
				$judge2="PASS:($ty1[0])-($ty2[0]):$hp1:$hp2:$pa1,$pb1,$pc1:$pa2,$pb2,$pc2";
			}
		}else{
			if(($ty1[0] eq "1|0" or $ty1[0] eq "0|1") and ($ty2[0] eq "1|0" or $ty2[0] eq "0|1")){
				$judge2="PASS:($ty1[0])-($ty2[0]):$hp1:$hp2:$pa1,$pb1,$pc1:$pa2,$pb2,$pc2";
			}else{
				$judge2="FAILED:($ty1[0])-($ty2[0]):$hp1:$hp2:$pa1,$pb1,$pc1:$pa2,$pb2,$pc2";
			}
		}
	}
	
	PHASE:
	$out->[1]=$judge2;
	
	#mapQ & PE
	my ($judge3,$judge4)=&judgemap($ref,$hfile,$index,$header);
	$out->[2]=$judge3;
	$out->[3]=$judge4;
	$out->[4]=$ex_line;

	# print "$r_bin\t$ref->[0]\t$ref->[1]\t$ref->[2]\t$ref->[3]\t$ref->[4]\t$ref->[5]\n$p12\t$p13\t$p14\t$p23\t$p24\t$p34\t$p_two\t$p_ex\n";
	return;
}


sub judgemap{
	my ($sv,$hfile,$index,$header)=@_;
	my $result;
	my @map=(0,0,0,0);
	my @pair=([()],[()],[()],[()]);
	
	&get_mapq_from_bam($sv->[2],$sv->[3],$sv->[4],$sv->[5],\@map,\@pair,$hfile,$index,$header);
	
	my $type=$sv->[6];
	my @count;
	$count[0]= scalar @{$pair[0]};
	$count[1]= scalar @{$pair[1]};
	$count[2]= scalar @{$pair[2]};
	$count[3]= scalar @{$pair[3]};
	my $count_sum=$count[0]+$count[1]+$count[2]+$count[3];
	my $pe_line;
	if($count_sum > 0){
		if($sv->[2] eq $sv->[4]){
			if($type eq "RL"){
				my $ratio=$count[0]/$count_sum;
				if($count[0] > 3 and $ratio >=0.75){
					$pe_line="PASS";
				}elsif($ratio < 0.75 and $count_sum > 3){
					$pe_line="FAILED";
				}else{
					$pe_line="NULL";
				}
			}elsif($type eq "LR"){
				my $ratio=$count[1]/$count_sum;
				if($count[1] > 3 and $ratio >=0.75){
					$pe_line="PASS";
				}elsif($ratio < 0.75 and $count_sum > 3){
					$pe_line="FAILED";
				}else{
					$pe_line="NULL";
				}
			}elsif($type eq "LL"){
				my $ratio=$count[2]+$count[3]/$count_sum;
				if($count[2] > 3 and $ratio >=0.75){
					$pe_line="PASS";
				}elsif($ratio < 0.75 and $count_sum > 3){
					$pe_line="FAILED";
				}else{
					$pe_line="NULL";
				}
			}elsif($type eq "RR"){
				my $ratio=$count[2]+$count[3]/$count_sum;
				if($count[3] > 3 and $ratio >=0.75){
					$pe_line="PASS";
				}elsif($ratio < 0.75 and $count_sum > 3){
					$pe_line="FAILED";
				}else{
					$pe_line="NULL";
				}
			}
		
		}else{
			if($type eq "RL"){
				my $ratio=($count[0]+$count[1])/$count_sum;
				if($count[0] > 3 and $ratio >=0.75){
					$pe_line="PASS";
				}elsif($ratio < 0.75 and $count_sum > 3){
					$pe_line="FAILED";
				}else{
					$pe_line="NULL";
				}
			}elsif($type eq "LR"){
				my $ratio=($count[0]+$count[1])/$count_sum;
				if($count[1] > 3 and $ratio >=0.75){
					$pe_line="PASS";
				}elsif($ratio < 0.75 and $count_sum > 3){
					$pe_line="FAILED";
				}else{
					$pe_line="NULL";
				}
			}elsif($type eq "LL"){
				my $ratio=($count[2]+$count[3])/$count_sum;
				if($count[2] > 3 and $ratio >=0.75){
					$pe_line="PASS";
				}elsif($ratio < 0.75 and $count_sum > 3){
					$pe_line="FAILED";
				}else{
					$pe_line="NULL";
				}
			}elsif($type eq "RR"){
				my $ratio=($count[2]+$count[3])/$count_sum;
				if($count[3] > 3 and $ratio >=0.75){
					$pe_line="PASS";
				}elsif($ratio < 0.75 and $count_sum > 3){
					$pe_line="FAILED";
				}else{
					$pe_line="NULL";
				}
			}
		}
	}else{
		$pe_line="NULL";
	}
	
	if($pe_line eq "NULL"){
		$pe_line=$pe_line.":$count[0],$count[1],$count[2],$count[3]:NA";
	}elsif($pe_line eq "FAILED"){
		$pe_line=$pe_line.":$count[0],$count[1],$count[2],$count[3]:NA";
	}else{
		if($type eq "RL"){
			my (@posa,@posb);
			foreach my $ref(@{$pair[0]}){
				push @posa,$ref->[0];
				push @posb,$ref->[1];
			}
			
			@posa=sort {$b <=> $a } @posa;
			@posb=sort {$a <=> $b } @posb;
			$pe_line=$pe_line.":$count[0],$count[1],$count[2],$count[3]:$posa[0]-$posb[0]";	
		}elsif($type eq "LR"){
			my (@posa,@posb);
			foreach my $ref(@{$pair[1]}){
				push @posa,$ref->[0];
				push @posb,$ref->[1];
			}
			
			@posa=sort {$a <=> $b } @posa;
			@posb=sort {$b <=> $a } @posb;
			$pe_line=$pe_line.":$count[0],$count[1],$count[2],$count[3]:$posa[0]-$posb[0]";	
		}elsif($type eq "LL"){
			my (@posa,@posb);
			foreach my $ref(@{$pair[2]}){
				push @posa,$ref->[0];
				push @posb,$ref->[1];
			}
			
			@posa=sort {$a <=> $b } @posa;
			@posb=sort {$a <=> $b } @posb;
			$pe_line=$pe_line.":$count[0],$count[1],$count[2],$count[3]:$posa[0]-$posb[0]";	
		}elsif($type eq "RR"){
			my (@posa,@posb);
			foreach my $ref(@{$pair[3]}){
				push @posa,$ref->[0];
				push @posb,$ref->[1];
			}
			
			@posa=sort {$b <=> $a } @posa;
			@posb=sort {$b <=> $a } @posb;
			$pe_line=$pe_line.":$count[0],$count[1],$count[2],$count[3]:$posa[0]-$posb[0]";	
		}
	}
		
	if($map[0] ==0 or $map[2] ==0){
		$result = "FAILED:NA-NA";
		goto PASS;
	}
	my $r1=sprintf("%.4f",$map[1]/$map[0]);
	my $r2=sprintf("%.4f",$map[3]/$map[2]);
	
	if($r1 < 0.5 and $r2 < 0.5){
		$result = "FAILED:$r1-$r2";
	}elsif($r1 < 0.5 or $r2 < 0.5){
		$result = "SUSPECT:$r1-$r2";
	}else{
		$result = "PASS:$r1-$r2";
	}
PASS:
	return ($result,$pe_line);
}


sub get_matrix{
	my ($chr1,$pos1,$chr2,$pos2,$extbin,$size,$ref1,$ref2,$hfile,$index,$header)=@_;
	
	my (@bar1,@bar2);
	my (@p1,@p2,@p3,@p4);
	my ($c1,$c2,$c3,$c4)=(0)x4;
	my %tmphash;
	#1
	my $start1=$pos1-$size/2;
	my $end1=$pos1+$size/2;
	$start1=0 if $start1 < 0;
	$end1=$reflen{$chr1} if $end1 > $reflen{$chr1};
	
	my $count=0;
	while($count < $extbin){
		my $new_start=$end1+$size*$count;
		my $new_end=$new_start+$size;
		last if $new_end > $reflen{$chr1};
		push @p1,$new_end;
		$count++;
	}
	
	$c1=$count;
	$count=0;
	while($count < $extbin){
		my $new_end=$start1-$size*$count;
		my $new_start=$new_end-$size;
		last if $new_start < 0; 
		unshift @p2,$new_start;
		$count++;
	}
	
	$c2=$count;
	$count=0;
	
	#2
	my $start2=$pos2-$size/2;
	my $end2=$pos2+$size/2;
	$start2=0 if $start2 < 0;
	$end2=$reflen{$chr2} if $end2 > $reflen{$chr2};
	
	while($count < $extbin){
		my $new_start=$end2+$size*$count;
		my $new_end=$new_start+$size;
		last if $new_end > $reflen{$chr2};
		push @p3,$new_end;
		$count++;
	}
	
	$c3=$count;
	$count=0;
	
	while($count < $extbin){
		my $new_end=$start2-$size*$count;
		my $new_start=$new_end-$size;
		last if $new_start < 0; 
		unshift @p4,$new_start;
		$count++;
	}
	
	$c4=$count;
	
	$count=$c1;
	$count=$c2 if $c2 < $count;
	$count=$c3 if $c3 < $count;
	$count=$c4 if $c4 < $count;
	
	return $count if $count <5;
	
	##
	my $n_cut;
	$n_cut=$c1-$count;
	for(my $i=0;$i<$n_cut;$i++){
		pop @p1;
	}
	
	$n_cut=$c2-$count;
	for(my $i=0;$i<$n_cut;$i++){
		shift @p2;
	}
	
	$n_cut=$c3-$count;
	for(my $i=0;$i<$n_cut;$i++){
		pop @p3;
	}
	
	$n_cut=$c4-$count;
	for(my $i=0;$i<$n_cut;$i++){
		shift @p4;
	}
	
	my $final_start1=$p2[0];
	my $final_end1=$p1[-1];
	
	my $final_start2=$p4[0];
	my $final_end2=$p3[-1];
	
	my @raw;
	&get_bar_from_bam($chr1,$final_start1,$final_end1,\@bar1,$hfile,$index,$header);
	&get_bar_from_bam($chr2,$final_start2,$final_end2,\@bar2,$hfile,$index,$header);
	
	@{$ref1}=@bar1;
	@{$ref2}=@bar2;
	
	return $count;
}

sub get_bar_from_bam{
	my ($chr,$start,$end,$ref,$hfile,$index,$header)=@_;
	
	for(my $i=$start;$i<$end;$i+=$size){
		my $index=int(($i-$start)/$size);
		%{$ref->[$index]}=();
	}
	
	my $callback = sub {
		my ($alignment,$data) = @_;
		# my $flag=$a->flag;
		#return if $flag & 0x400;
		my $name=$alignment->qname;
		my $barcode=(split(/#/,$name))[-1];
		return if $barcode eq "0_0_0";
		my $pos=($alignment->pos)+1;
		return if $pos < $start;
		return if $pos > $end;
		
		my $index=int(($pos-$start)/$size);
		$ref->[$index]->{$barcode}=1;
	};
	$index->fetch($hfile,$header->parse_region("$chr:$start-$end"),$callback);
	
	return;
}


# sub get_bar_from_bam1{
	# my ($chr,$start,$end,$ref)=@_;

	# for(my $i=$start;$i<$end;$i+=$size){
		# my $index=int(($i-$start)/$size);
		# %{$ref->[$index]}=();
	# }
	
	# my $fh=IO::Pipe->new();
	# $fh->reader("samtools view -F 1028 -@ 2 $bamfile $chr:$start-$end");
	# # open my $fh,"samtools view -F 1028 -@ 2 $chr.bam $chr:$start-$end|" or die $!;
	# while ( ! ($fh->eof) ) {
		# my $line;
		# defined($line = $fh->getline) or die "readline failed at $chr:$start-$end: $!";
		# chomp $line;
		# my @t=split(/\t/,$line);
		# my $name=$t[0];
		# my $barcode=(split(/#/,$name))[-1];
		# next if $barcode eq "0_0_0";
		# my $pos=$t[3];
		# print STDERR "$line\n" unless (defined $pos);
		# next if $pos < $start;
		# next if $pos > $end;
		# my $index=int(($pos-$start)/$size);
		# $ref->[$index]->{$barcode}=1;
	# }
	# $fh->close;
	
	# return;
# }

sub get_mapq_from_bam{
	
	my ($chr1,$pos1,$chr2,$pos2,$map,$pe,$hfile,$index,$header)=@_;
	my $N_search;
	
	if($chr1 eq $chr2){
		$N_search=int((($pos2-$pos1)/$size/4)+0.5);
		$N_search=$Nmerge if $N_search > $Nmerge;
	}else{
		$N_search=$Nmerge;
	}
	$N_search =1 if $N_search < 1;
	
	my $start1=$pos1-$size*$N_search;
	my $end1=$pos1+$size*$N_search;
	$start1=0 if $start1 <0;
	
	my $start2=$pos2-$size*$N_search;
	my $end2=$pos2+$size*$N_search;
	$start2=0 if $start1 <0;
	
	my $callback1 = sub {
		my ($alignment,$data) = @_;
		my $flag=$alignment->flag;
		return if $flag & 0x400;
		my $pos=($alignment->pos)+1;
		return if $pos < $start1;
		return if $pos > $end1;
		
		my $q=$alignment->qual;
		$map->[0]++;
		if($q > 10){
			$map->[1]++;
		}
		
		return if $flag & 0x04;
		return if $flag & 0x08;
		
		my $mtid=$alignment->mtid;
		my $mchr=$header->target_name->[$mtid];
		return if $mchr ne $chr2;
		
		my $mpos=($alignment->mpos)+1;
		return if $mpos < $start2;
		return if $mpos > $end2;
		
		my $type;
		if($flag & 0x10){
			if($flag & 0x20){
				$type="--";
			}else{
				$type="-+";
			}
		}else{
			if($flag & 0x20){
				$type="+-";
			}else{
				$type="++";
			}
		}
		
		if($type eq "+-"){
			push @{$pe->[0]},[($pos,$mpos)];
		}
		
		if($type eq "-+"){
			push @{$pe->[1]},[($pos,$mpos)];
		}
		
		if($type eq "--"){
			push @{$pe->[2]},[($pos,$mpos)];
		}
		
		if($type eq "++"){
			push @{$pe->[0]},[($pos,$mpos)];
		}	
	};
	$index->fetch($hfile,$header->parse_region("$chr1:$start1-$end1"),$callback1);
	
	my $callback2 = sub {
		my ($alignment,$data) = @_;
		my $flag=$alignment->flag;
		return if $flag & 0x400;
		my $pos=($alignment->pos)+1;
		return if $pos < $start2;
		return if $pos > $end2;
		
		my $q=$alignment->qual;
		$map->[2]++;
		if($q > 10){
			$map->[3]++;
		}
		# usleep(10);
	};
	$index->fetch($hfile,$header->parse_region("$chr2:$start2-$end2"),$callback2);
	return;
}

# sub get_mapq_from_bam1{
	
	# my ($chr1,$pos1,$chr2,$pos2,$map,$pe,$svid)=@_;
	# my $start1=$pos1-$size*$Nmerge;
	# my $end1=$pos1+$size*$Nmerge;
	# $start1=0 if $start1 <0;
	
	# my $start2=$pos2-$size*$Nmerge;
	# my $end2=$pos2+$size*$Nmerge;
	# $start2=0 if $start2 <0;
	
	# my $fh1=IO::Pipe->new();
	# $fh1->reader("samtools view -F 1028 -@ 2 $bamfile $chr1:$start1-$end1");
	# # open my $fh1,"samtools view -F 1028 -@ 2 $chr1.bam $chr1:$start1-$end1|" or die $!;
	# print "test1\n";
	# while(! ($fh1->eof)){
		# my $line;
		# defined($line = $fh1->getline) or die "readline failed at $chr1:$start1-$end1: $!";
		# chomp $line;
		# my @t=split(/\t/,$line);
		# my $pos=$t[3];
		# print STDERR "$_\n" unless (defined $pos);
		# next if $pos < $start1;
		# next if $pos > $end1;
		# my $q=$t[4];
		# $map->[0]++;
		# if($q > 10){
			# $map->[1]++;
		# }
		
		# my $flag=$t[1];
		# # next if $flag & 0x04;
		# next if $flag & 0x08;
		# my $mchr=$t[6];
		# if($mchr eq "="){
			# $mchr=$chr1;
		# }
		# next if $mchr ne $chr2;
		# my $mpos=$t[7];
		# next if $mpos < $start2;
		# next if $mpos > $end2;
		
		# my $type;
		# if($flag & 0x10){
			# if($flag & 0x20){
				# $type="--";
			# }else{
				# $type="-+";
			# }
		# }else{
			# if($flag & 0x20){
				# $type="+-";
			# }else{
				# $type="++";
			# }
		# }
		
		# if($type eq "+-"){
			# push @{$pe->[0]},[($pos,$mpos)];
		# }
		
		# if($type eq "-+"){
			# push @{$pe->[1]},[($pos,$mpos)];
		# }
		
		# if($type eq "--"){
			# push @{$pe->[2]},[($pos,$mpos)];
		# }
		
		# if($type eq "++"){
			# push @{$pe->[0]},[($pos,$mpos)];
		# }
	# }
	# $fh1->close;
	# print "S3-1-1\t$svid\n";
	
	# my $fh2=IO::Pipe->new();
	# $fh2->reader("samtools view -F 1028 -@ 2 $bamfile $chr2:$start2-$end2");
	# # open my $fh2,"samtools view -F 1028 -@ 2 $chr2.bam $chr2:$start2-$end2|" or die $!;
	# while(! ($fh2->eof)){
		# my $line;
		# my $t0 = [gettimeofday];
		# defined($line = $fh2->getline) or die "readline failed at $chr2:$start2-$end2: $!";
		# chomp $line;
		# my @t=split(/\t/,$line);
		# my $pos=$t[3];
		# next if $pos < $start2;
		# next if $pos > $end2;
		# my $q=$t[4];
		# $map->[2]++;
		# if($q > 10){
			# $map->[3]++;
		# }
		# my $elapsed = tv_interval ($t0); 
		# print "$elapsed\n";
	# }
	# $fh2->close;
	# print "S3-1-2\t$svid\n";
	
	# return;
# }


sub get_type{
	my ($id,$ref)=@_;
	foreach my $i(@{$ref}){
		my @tmp=split(/:/,$i);
		
		#jump the exists
		next if exists $link{$id}{$tmp[0]};
		
		#assign sv id
		$svid++;
		#cache the id to sv
		$link{$id}{$tmp[0]}=$svid;
		$link{$tmp[0]}{$id}=$svid;
		
		my (@info0,@info1);
		if(exists $clusterinfo{"$id-$tmp[0]"}){
			@info0=@{$clusterinfo{"$id-$tmp[0]"}}[0..2];
			@info1=@{$clusterinfo{"$id-$tmp[0]"}}[3..5];
		}elsif(exists $clusterinfo{"$tmp[0]-$id"}){
			@info1=@{$clusterinfo{"$tmp[0]-$id"}}[0..2];
			@info0=@{$clusterinfo{"$tmp[0]-$id"}}[3..5];
		}else{
			die "Wrong linked all file!!!\n";
		}
		
		##creat the SV
		my @local_sv;
		if($info0[0] eq $info1[0]){
			if($info0[1] < $info1[1]){
				my $chr1=$info0[0];
				my $pos1=$info0[1];
				my $chr2=$info1[0];
				my $pos2=$info1[1];
				my $type=$info0[2].$info1[2];
				@local_sv=($id,$tmp[0],$chr1,$pos1,$chr2,$pos2,$type,$tmp[1]);
			}else{
				my $chr1=$info1[0];
				my $pos1=$info1[1];
				my $chr2=$info0[0];
				my $pos2=$info0[1];
				my $type=$info1[2].$info0[2];
				@local_sv=($tmp[0],$id,$chr1,$pos1,$chr2,$pos2,$type,$tmp[1]);
			}
		}else{
			my $nchr1=$info0[0];
			my $nchr2=$info1[0];
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
				my $chr1=$info0[0];
				my $pos1=$info0[1];
				my $chr2=$info1[0];
				my $pos2=$info1[1];
				my $type=$info0[2].$info1[2];
				@local_sv=($id,$tmp[0],$chr1,$pos1,$chr2,$pos2,$type,$tmp[1]);
			
			}else{
				my $chr1=$info1[0];
				my $pos1=$info1[1];
				my $chr2=$info0[0];
				my $pos2=$info0[1];
				my $type=$info1[2].$info0[2];
				@local_sv=($tmp[0],$id,$chr1,$pos1,$chr2,$pos2,$type,$tmp[1]);
			
			}	
		}	
		$svinfo{$svid}=[@local_sv];
	}
	return;
}

sub search_extend{
	my $ref=$_[0];
	
	my $max="NULL";
	my $index;
	for(my $i=0;$i<$#{$ref};$i++){
		my $before=0;
		for(my $j=0;$j<=$i;$j++){
			$before+=$ref->[$j];
		}
		
		my $after=0;
		for(my $j=$i+1;$j<=$#{$ref};$j++){
			$after+=$ref->[$j];
		}
		
		my $value=$before-$after;
		if($max eq "NULL"){
			$max=$value;
			$index=$i;
		}else{
			if($value >= $max){
				$max=$value;
				$index=$i;
			}
		}
	}
	
	return $index+1;
}
###############################################################






