#!usr/bin/perl -w
use strict;
use Getopt::Long;
use PerlIO::gzip;

############################################################################################################
#### This program is to cluster the discrodant reads(paired reads) that support the breakpoint			####
#### Author: zeng xi																					####
#### Last update: 2020-12																				####
############################################################################################################

my $usage="perl $0 -i <se_se file> -h <human_un file> -l <gap length> -o <output>";
my ($in,$len,$human,$out);
GetOptions (
	'i=s'=>\$in,
	'l=i'=>\$len,
	'o=s'=>\$out,
	'h=s'=>\$human,
);
die $usage if !$in|!$len|!$out|!$human;
my %breakpoint;
if ($in=~/.gz\b/){
	open IN,"gzip -dc $in|" or die $!;
}else{
	open IN,"$in" or die $!;
}
while(<IN>){
	chomp;
	my ($chr,$start)=(split /\s+/)[1,2];
	$breakpoint{$chr}{$start}++;
	<IN>;

}
close IN;

my %breakhuman;
open HU,"gzip -cd $human|" or die "can't open $human\n";
while (<HU>){
	chomp;
	my ($chr,$start)=(split /\s+/)[0,1];
	for my $pos (keys %{$breakpoint{$chr}}){
		if($start<=$pos+$len && $start>=$pos-$len){
			$breakhuman{$chr}{$pos}++;
			last;
		}
	}
	<HU>;
}
close HU;
	
my ($num,$human_num,$start)=(0,0,0);
my ($old_chr,$old_postion,@array);
open OUT,">$out" or die $!;
print OUT "chr\tbreakpoint\treads num\thunman num\ttotal num \n";
for my $chr (keys %breakpoint){
	for my $postion(sort {$a<=>$b}keys %{$breakpoint{$chr}}){
		if ($postion-$start>$len || $old_chr ne $chr){
			$old_postion= &medain(@array) if @array!=0;
			print OUT "$old_chr\t$old_postion\t$num\t$human_num\t",$num+$human_num,"\n"  unless $num==0;
			$start=$postion;
			$num=0;
			$human_num=0;
			@array=();
			$breakpoint{$chr}{$postion}=0 if !exists $breakpoint{$chr}{$postion};
			$num=$breakpoint{$chr}{$postion}+$num;
			$breakhuman{$chr}{$postion}=0 if !exists $breakhuman{$chr}{$postion};
			$human_num=$breakhuman{$chr}{$postion}+$human_num;
			push (@array,$postion);
			$old_chr=$chr;
		}else{
			$breakpoint{$chr}{$postion}=0 if !exists $breakpoint{$chr}{$postion};
			$num=$breakpoint{$chr}{$postion}+$num;
			$breakhuman{$chr}{$postion}=0 if !exists $breakhuman{$chr}{$postion};
			$human_num=$breakhuman{$chr}{$postion}+$human_num;
			push (@array,$postion);
			$old_chr=$chr;
		}
	}
}
close OUT;
sub medain{
	my @data=@_;
	return($data[int(@data/2)]);
}
			
