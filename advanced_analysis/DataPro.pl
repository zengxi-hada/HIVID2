#!/usr/bin/perl

use strict;
use warnings;

my $data_file = shift;

die "perl $0 <data_file>\n" unless ($data_file);


my %h;
open DATA,$data_file or die $!;

while(<DATA>){
	chomp;
	next if /gene/;
	next if /^\s+/;
	my @a = split;
	my $lib = $a[0];
	$h{$lib}++;
}
print "GENE_NAME\t\t\t\tNumber\n";
for my $lib(keys %h){
	printf "$lib\t\t\t\t$h{$lib}\n";
}
