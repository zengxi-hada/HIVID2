#!/usr/bin/perl

use strict;
use warnings;

my $infile = shift;
my $out = shift;
my $threshold = shift;
my $filter = shift;

die "perl $0 <mutual sam file> <initialized_bk> <threshold of alignment length in bwasw> <filter nouniq align or not(y|n)>\n" unless ($infile && $out && $threshold && $filter);

#######################################################################################################################
######## this program is to assign right_support or left_support to each support read of virus integration  ###########
######## Author: zeng xi																					###########
######## Last update: 2020-12																				###########
#######################################################################################################################

open IN, $infile or die $!;
open OUT, ">$out" or die $!;

print OUT "ref\tpos\t\treads id\n";

while(<IN>){
	next if /reads_id/; 
	chomp;
	my @a = split;
	my $cigar = $a[2];
	my $ref = $a[1];
	my $pos = $a[3];
	my $id = $a[0];
	if($filter eq "y"){
		next if ($id =~ /repeat_reads/);
	}
	my $index = 0;
	my $precise_pos = 0;
	my $m;
	my $orientation;
	while($cigar =~ /(\d+)M/g){
        $m += $1;
    }
    if($m < $threshold){next;}
	my $num; my $num2; my $num1;my $num_tmp;
	if(($num) = $cigar=~ /^(\d+)S.*\d+M$/){
		if($num >= 5){$orientation = "right_support";}							# the read is located on the downstream of the virus integration
	}elsif(($num_tmp, $num) = $cigar=~ /^\d+M(.*\D)*(\d+)S$/){
		if($num >= 5){$orientation = "left_support";}							# the read is located on the downstream of the virus integration
	}elsif((($num)= $cigar=~/^.*(\d+)I.*$/) && ($num >= 15)){					# the read contains a virus integration in the middle of itself ?
		$orientation = "left_right";
	}elsif(($num1, $num_tmp, $num2) = $cigar=~ /^(\d+)S(.*\D)*(\d+)S$/){
		if($num1 >= 5 && $num2 >= 5){
			$orientation = "left_right";										# the read contain virus integration on both side of itself
		}elsif($num1 >= 5 && $num2 < 5){
			$orientation = "right_support";
		}elsif($num1 < 5 && $num2 >= 5){
			$orientation = "left_support";
		} 
	}	
	if(!$orientation){$orientation = "unrecognised";}
	while($cigar =~ /((\d+)[MSID])/g){
		my $len = $2;
		if($1 =~ /S/){
			next if $len < 14;													# updated at 2019-07-02
			$precise_pos = $index + $pos;
			print OUT "$ref\t$precise_pos\t$orientation\t$cigar\t$id\n";
			next;
		}elsif($1 =~ /I/){
			next if $len < 14;													# updated at 2019-07-03
			$precise_pos = $index + $pos;
			print OUT "$ref\t$precise_pos\t$orientation\t$cigar\t$id\n";
			next;
		}
		$index += $len;
	}
}

close IN;
