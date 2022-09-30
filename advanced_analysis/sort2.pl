#!usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $file =shift;
chomp (my $out_file=`pwd`);
print $out_file."\n";
$out_file.="/all_bk_sort.anno";
print $out_file."\n"; 
open IN,$file or die $!;
open OUT,">$out_file" or die $!;
print OUT "sample_name\tchr\tpos\tgene\tregion\n";
$|=1;
while(<IN>)
{
   chomp;
   next if /ref/;
   my @a=split;
   my $name=$a[0];
   my $chr=$a[1];
   my $pos=$a[2];
   (my $region)=$a[3]=~m/(intergenic|intron|5-UTR|3-UTR|CDS)/;
   my @b;
   push @b,$a[3]=~m/(\w+);:/;
   my $i=4;
   while(exists($a[$i]))
   {
     # (my @gene)=$a[4]=~m/(\w+);:/;
     # print "@gene\n";  
      push @b,$a[$i]=~m/(\w+);:/;
      $i++;
      #push @b,$a[4]=~m/*;:*/;
   }
   my $b;
   for my $i(@b)
   {
      $b.=$i.",";
   }
   print "$b\n";
   print OUT "$name\t$chr\t$pos\t$b\t$region\n";
}
