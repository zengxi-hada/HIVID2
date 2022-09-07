#!/usr/bin/perl -lan

use strict;
use warnings;

my $spos = $F[1];
my $epos = $F[2];
my $chr = $F[0];
my $cluster_value = $F[3];

my ($chr_num) = $chr=~/chr(\w+)$/;

#my $spos = $pos*1000;
#my $epos = $spos + 999; 

print "hs$chr_num\t$spos\t$epos\t$cluster_value"; 
