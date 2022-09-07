#!/usr/bin/perl

use strict;
use warnings;

my $input = shift;


die "perl $0 <input breakpoint file>\n" unless ($input);

open IN, $input or die $!;

my %h;
my @isample;
my @ipos;

while(<IN>){
	next if /ref/; 
	chomp;
	my @a = split;
	my $sample = $a[0];
	my $chr = $a[1];
	my $pos = $a[2];
	my $support = $a[3];
	my $combination = "$sample\t$chr\t$pos";
	$h{$combination} = $support;
	push @isample, $sample;
	push @ipos, "$chr\t$pos";
}
close IN;

my %count;
@isample = grep {++$count{$_} < 2} @isample;


my $support;
my $sum;
my $average_sup;
my $flag = 0;

=h
for my $i2(@isample){
	$sum = 0;
	if($flag == 0){
		for my $i1(@ipos){
			my $chr = (split /\s+/, $i1)[0];
			my $pos = (split /\s+/, $i1)[1];
			print "$chr\_$pos\t";
		} 
		$flag = 1;
		print "\n";
	}

	for my $i1(@ipos){
		my $combination = "$i2\t$i1";
		if(exists $h{$combination}){
			printf "%.3f\t", $h{$combination};
		}else{
			print 0, "\t";
		}
	}
	print "\n";
}
=cut

for my $i2(@ipos){
	$sum = 0;

=h
	if($flag == 0){
		for my $i1(@isample){
			print "$i1\t";
		}
		$flag = 1;
		print "\n";
	}
=cut
	
	my $chr = (split /\s+/, $i2)[0];
	my $pos = (split /\s+/, $i2)[1];
	print "$chr\_$pos\t";
	my $count = 0;
	for my $i1(@isample){
		my $combination = "$i1\t$i2";
		if(exists $h{$combination}){
			printf "%.3f\t", $h{$combination};
			$count += $h{$combination};
#			printf  "$h{$combination}\t";
		}else{
			print 0, "\t";
		}
	}
	print "$count";
	print "\n";
}

