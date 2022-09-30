#!usr/bin/perl

use strict;
use warnings;

my $file =shift;

die "perl $0 <file>\n" unless $file;
my %m=(
	   '1'=>'chr1',
       '2'=>'chr10', 
       '3'=>'chr11',
       '4'=>'chr12', 
       '5'=>'chr13', 
       '6'=>'chr14', 
       '7'=>'chr15', 
       '8'=>'chr16',
       '9'=>'chr17', 
       '10'=>'chr18', 
       '11'=>'chr19',
       '12'=>'chr2', 
       '13'=>'chr20', 
       '14'=>'chr21', 
       '15'=>'chr22', 
       '16'=>'chr3', 
       '17'=>'chr4', 
       '18'=>'chr5', 
       '19'=>'chr6', 
       '20'=>'chr7', 
       '21'=>'chr8',
       '22'=>'chr9'
);

open FILE,$file or die $!;

my %h;
while(<FILE>){
	chomp;
	next if /chr/;
	my @a =split;
	#if($a[1]!=0) {$h{'chr1'}{'count'}++;}
	my $i;
	for($i=1;$i<=44;$i+=2)
	{
		if($a[$i]!=0){
			$h{$i/2+0.5}{'num'}+=$a[$i];
		}
		$h{$i/2+0.5}{'ratio'}+=$a[$i+1];
	}
}
#my $except=1.10111;
#my $exceptRatio=0.020390926;
my $j;
my $num1;
for($j=1;$j<=44;$j+=2){
    my $num = $h{$j/2+0.5}{'num'};
    #my $ratio = $h{$j}{'ratio'};
	$num1+=$num;
	#$exceptRatio+=$ratio;    
}

my $i;
my $exceptRatio;
if($num1!=0)
{
   $exceptRatio=$h{1}{'ratio'}/$num1;
}
else {$exceptRatio=0;}
for($i=1;$i<=44;$i+=2){
	my $r=$i/2+0.5;
	my $num;
	if(exists($h{$r}{'num'}) && $h{$r}{'num'} ne ""){
		$num = $h{$r}{'num'};
	}
	else {$num =0;}
	my $except = $h{$r}{'ratio'};
	my $chrnum = $m{$r};
	my $ratio;
 	if($num1!=0)
	{$ratio=$num/$num1;}	
	else {$ratio=0;}
	print "$chrnum\t $num \t $ratio  \t$except\t$exceptRatio\n";
}
close FILE;
