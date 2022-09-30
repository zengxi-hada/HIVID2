#!usr/bin/perl

use strict;
use warnings;
use File::Basename;

#my $add_file=shift;
my @a=`ls *.final`;
my $out=`pwd`;
chomp($out);
$out.="/break_all.txt";
print $out."\n";
open OUT,">$out" or die $!;
print OUT "name\tref\tpos\n";
foreach my $i(@a)
{
   my $name="";
   open IN,$i or die $!;
   ($name)=$i=~m/(ERS\d+)*/;
   print "$name\n";
  # if($i=~m/(ERS\d+)/) {print "111\n";$name=$1;}
   while(<IN>)
   {
      next if /ref/;
      chomp;
     # if($i=~m/(ERS*)_/) {$name=$1;}
      my @b=split;
      my $chr=$b[0];
      my $pos=$b[1];
      print OUT "$name\t$chr\t$pos\n";
   }
   close IN;	
}
#close OUT;
