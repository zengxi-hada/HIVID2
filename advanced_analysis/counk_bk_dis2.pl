#!/sur/bin/perl -w
use strict;

my $alu = shift;
my $bk = shift;
my $out = shift;

die "perl $0 <Alu.txt> <human.Tumor.bk.gte10> <bk_Alu.dis>" unless($alu and $bk and $out);

open ALU, $alu or die $1;
open BK, $bk or die $!;
open OUT, ">$out" or die $!;

my %hash;
while(<ALU>){
	chomp;
	my @tmp = split;
	$hash{$tmp[0]}{$tmp[1]} = $tmp[2];
}
close ALU;

while(<BK>){
	chomp;
	my @tmp =split;
	my $pos = $tmp[2];
	my $dis_tmp = 0;
	my $dis = "nan";
	my $spos;
	my $epos;
	for my $start (sort {$a<=>$b} keys %{$hash{$tmp[1]}}){
		my $end = $hash{$tmp[1]}{$start};
		if($dis_tmp == 0){
			if ($start > $pos){
				$dis_tmp = $start - $pos;
				$dis = $dis_tmp;
				$spos = $start;	$epos = $end;
#print "AA\t$dis_tmp\n";
				last;
			}elsif($start <= $pos && $end >= $pos){
				$dis = 0;
				$spos = $start; $epos = $end;
#print "BB\t$dis_tmp\n";
				last;
			}else{
				$dis_tmp = $pos - $end;
#print "BB1\t$dis_tmp\n";
			}
		}else{
			if($start > $pos){
				my $dis_tmp2;
				$dis_tmp2 = $start - $pos;
				if($dis_tmp2 < $dis_tmp ){
					$dis_tmp = $dis_tmp2;
				}
				$dis = $dis_tmp;
				$spos = $start; $epos = $end;
#print "CC\t$dis_tmp\n";
				last;
			}elsif($start <= $pos && $end >= $pos){
				$dis = 0;
				$spos = $start; $epos = $end;
#print "DD\t$dis_tmp\n";
				last;
			}else{
				$dis_tmp = $pos - $end;
				$dis = $dis_tmp;
				$spos = $start; $epos = $end;
#print "DD1\t$dis_tmp\n";
			}
		}
#		if($dis eq "nan"){
#			$dis = $dis_tmp;
#$spos = $start; $epos = $end;
#print "EE\t$dis_tmp\n";
#		}
	}
#print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$dis\t$spos\t$epos\n";
	if($dis ne "nan"){
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$dis\t$spos\t$epos\n";
	}
}
