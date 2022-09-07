#!/usr/bin/perl

use strict;
use warnings;

my $infile = shift;
my $out = shift;
my $threshold = shift;
my $filter = shift;

die "perl $0 <mutual sam file> <initialized_bk> <threshold of alignment length in bwasw> <filter nonuniq reads or not(y|n)>\n" unless ($infile && $out && $threshold && $filter);

########################################################################################################
# this program is to assign right_support or left_support to each support read of virus integration ####
# Author: zeng xi																					####
# Last update: 2020-12																				####
#### ###################################################################################################

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
	my $num; my $num2; my $num1; my $num_tmp;
    if(($num) = $cigar=~ /^(\d+)S.*\d+M$/){
        if($num >= 5){$orientation = "right_support";}
    }elsif(($num_tmp, $num) = $cigar=~ /^\d+M(.*\D)*(\d+)S$/){
		if($num >= 5){$orientation = "left_support";}
	}elsif((($num) = $cigar=~/^.*(\d+)I.*$/) && ($num >= 5)){
		$orientation = "left_right";
    }elsif(($num1,$num_tmp, $num2) = $cigar=~ /^(\d+)S(.*\D)*(\d+)S$/){
        if($num1 >= 5 && $num2 >= 5){
            $orientation = "left_right";
        }elsif($num1 >= 5 && $num2 < 5){
            $orientation = "right_support";
        }elsif($num1 < 5 && $num2 >= 5){
            $orientation = "left_support";
        }

    }
    if(!$orientation){$orientation = "unrecognised";}

	my $flags = 0;
	my $flagi = 0;
	while($cigar =~ /((\d+)[MSID])/g){
		my $len = $2;
		my $type = $1;
		if($type =~ /S/){
			next if $len < 14;												# updated at 2019-07-02
			$precise_pos = $index + $pos;
			print OUT "$ref\t$precise_pos\t$orientation\t$cigar\t$id\n";
			$flags = 1;
			next;
		}elsif($type =~ /(\d+)I/){
			next if $len < 14;												# updated at 2019-07-03
			my $insert_len = $1;
			my $precise_pos1 = $index + $pos;
			print OUT "$ref\t$precise_pos1\t$orientation\t$cigar\t$id\n";
			$flagi = 1;
			next;
		}
		$index += $len;
	}
}

close IN;
