#!/usr/bin/perl

use strict;
use warnings;
#use PerlIO::gzip;
use Getopt::Long;
#use File::Basename;
#use Cwd qw/abs_path/;
#use FindBin qw($Bin $Script);
#use File::Basename qw(basename dirname);
#use Data::Dumper;
#use File::Path qw(make_path remove_tree);

=head1 usage

  perl $0
		-bk	<str>		break point file
		-l	<str>		chr length file
		-w	<int>		size of window [1000]
		-ol <int>		overlap of two neighbouring windows
		-o	<str>		out put file

=cut

my ($bk_file, $chr_length, $out, $overlap);
my $window = 1000;

GetOptions(
	'bk=s' => \$bk_file,
	'l=s' => \$chr_length,
	'w=i' => \$window,
	'ol=i' => \$overlap,
	'o=s' => \$out,
);

die `pod2text $0` unless ($bk_file && $chr_length && $out && $overlap);

if($overlap == 1){$overlap = 0;}

my (@chr_length, %h);

open LE, $chr_length or die $!;
while(<LE>){
   chomp;
   next if /Total/;
   my @a = split /\s+/, $_;
   my $tmp = "$a[0]\t$a[1]";
   push @chr_length, $tmp;
}
close LE;

for my $i(@chr_length){
   my @split = split /\t/, $i;
   my $chr = "$split[0]";
   my ($start, $end);
   for (my $ipos=0; $ipos<=$split[1]; $ipos += ($window-$overlap)){
	  $start = $ipos+1;
	  $end = $ipos + $window;
      my $ikey = "$start\t$end";
      $h{$ikey}{"count"} = 0;
   }
}

my %record;
open BK, $bk_file or die $!;
while(<BK>){
	chomp;
	next if /ref/;
	my @a = split;
	my $sample = $a[0];
	my $chr = $a[1];
	my $pos = $a[2];
	my ($ichr, $istart, $iend);
#	print "zengxi\t$pos\n";
	for my $i(keys %h){
		$istart = [split /\t/, $i] -> [0];
		$iend = [split /\t/, $i] -> [1];
#		print "istart#$istart\tiend#$iend\n";
##		$ichr = [split /\t/, $i] -> [0];
##		if($pos>=$istart && $pos<=$iend && ($ichr eq $chr)
#:q
		if($pos=~m/(\d+)\(\d+\)/) { $pos=$1;}
		if($pos>=$istart && $pos<=$iend && (not exists $record{"$istart\t$iend"}{"record"}{$sample})){
#			print "$sample\n";
			$h{"$istart\t$iend"}{"count"}++;
			$record{"$istart\t$iend"}{"record"}{$sample} = 1;
			last;
		}
	}
}
close BK;

open OUT, ">$out" or die $!;
for my $i(sort keys %h){
	my $print = join "\t", (split /\t/, $i), "\t", $h{$i}{"count"}, "\n";
	print OUT "HBV_C\t$print";
}
close OUT;
