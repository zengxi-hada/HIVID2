#!/usr/bin/perl -w
use strict;
use lib '/public/home/xzeng/bin/BGI_bin/basic_soft/localperl/lib/';
use Cwd qw/abs_path/;
use lib '/home/xzeng';
use lib '/public/home/xzeng/.cpan/build';
use Getopt::Long;
use File::Path;
#use talmud;

my $usage="
	
Options:
                        -i <STR>         The distance files input;
                        -o <STR>         The output dile [current dir];
#                       -w <INT>         the length of windows(all 100 windows);
                        -h|?             Help!

                Example:
                        perl $0 -i <tumor.bk.Alu.dis> -w <100> -o <tumor.dis.num.win>
";
my ($in,$win,$out,$help);
GetOptions(
        'i=s' => \$in,
        'w=i' => \$win,
        'o=s' => \$out,
        'h|?' => \$help
);


die "$usage\n" if ($help || !$in);
$out=abs_path($out);

open IN, $in or die $!;
open OUT, ">$out" or die $!;

my %hash;	my %h;
my $h0 = 0;
#for (0..99){
for (0..19.5*100/$win){
	my $start;
	my $end;
	$start = $_ * $win + 1;
	$end = $start + $win - 1;
	$hash{$start,$end} = 0;
	$h{$start} = $end;
#print "$start\_$end\n";
}
my $bk;
while(<IN>){
	$bk++;
	chomp;
	my @tmp = split;
	if($tmp[3] == 0){
		$h0++;
	}
	for my $sta (sort {$a<=>$b} keys %h){
#print "$sta\t$h{$sta}\t$tmp[3]\n";
		if($sta<=$tmp[3] && $h{$sta}>=$tmp[3]){
			$hash{$sta,$h{$sta}}++;
#print "$sta\t$h{$sta}\t$tmp[2]\t$hash{$sta,$h{$sta}}\n";
		}
	}
}
my $h0_tmp;
$h0_tmp = $h0/$bk*1000;
print OUT "0\t0\t$h0_tmp\n";
for my $st (sort {$a<=>$b} keys %h){
	my $bbb;
	$bbb = $hash{$st,$h{$st}}/$bk*1000;
#	print OUT "$st\t$h{$st}\t$hash{$st,$h{$st}}\n";
	print OUT "$st\t$h{$st}\t$bbb\n";
}
