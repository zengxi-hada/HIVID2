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
   my $tmp = "$a[0]_$a[1]";
   push @chr_length, $tmp;
}
close LE;

for my $i(@chr_length){
   my @split = split /_/, $i;
   my $chr = "$split[0]";
   my ($start, $end);
   for (my $ipos=0; $ipos<=$split[1]; $ipos += ($window-$overlap)){
	  $start = $ipos+1;
	  $end = $ipos + $window;
      my $ikey = "$chr\_$start\_$end";
      $h{$ikey}{"count"} = 0;
   }
}

open BK, $bk_file or die $!;
while(<BK>){
	chomp;
	next if /ref/;
	my @a = split;
	my $sample = $a[0];
	my $chr = $a[1];
	my $pos = $a[2];
#	my $chr_name;
#	($chr_name) = $chr=~ /HBV_(\w)\d/;
	my ($ichr, $istart, $iend);
	for my $i(keys %h){
		$istart = [split /_/, $i] -> [1];
		$iend = [split /_/, $i] -> [2];
		$ichr = [split /_/, $i] -> [0];
#		print "$pos\t$istart\t$iend\n";
		if($pos>=$istart && $pos<$iend && ($ichr eq $chr && not exists $h{"$chr\_$istart\_$iend"}{"record"}{$sample})){
			$h{"$chr\_$istart\_$iend"}{"count"}++;
			$h{"$chr\_$istart\_$iend"}{"record"}{$sample} = 1;
			last;
		}
	}
}
close BK;


open OUT, ">$out" or die $!;
for my $i(sort keys %h){
	my $print = join "\t", (split /_/, $i), "\t", $h{$i}{"count"}, "\n";
	print OUT "$print";
}
close OUT;
