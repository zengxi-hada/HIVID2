#!/usr/bin/perl -lan

use strict;
use warnings;

my $pos = [split /_/, $F[0]]->[1];
my $chr = [split /_/, $F[0]]->[0];
my $heat_value = $F[-1];

my ($chr_num) = $chr=~/chr(\w+)$/;


print "hs$chr_num\t$pos\t$pos\t$heat_value";
 
