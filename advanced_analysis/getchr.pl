#!usr/bin/perl
use warnings;
use strict;


my %ratio = ('chr1' => 0.0080,
       'chr2' => 0.0082,
        'chr3' => 0.0101,
        'chr4' => 0.0105,
        'chr5' => 0.0109,
        'chr6' => 0.0116,
        'chr7' => 0.0125,
                'chr8' =>0.0137,
        'chr9' => 0.0144,
        'chr10' => 0.0149,
        'chr11' => 0.0148,
        'chr12' => 0.0150,
        'chr13' => 0.0175,
        'chr14' => 0.0186,
        'chr15' => 0.0196,
        'chr16' => 0.0222,
        'chr17' => 0.0240,
        'chr18' => 0.0250,
        'chr19' => 0.0338,
        'chr20' => 0.0312,
        'chr21' => 0.0425,
        'chr22' => 0.0392,
);
=cut
my $normal=shift;
my $tumor=shift;

die "<normal> <tumor>" unless ($normal && $tumor);
my (%h1,%h2);
open Nor,$normal or die $!;
open Tum,$tumor or die $!;

my ($count1,$count2);
while(<Nor>)
{
  chomp;
  next if /X/;
  next if /sample/;
  my @a=split;
  my $chr=$a[1];
  $h1{$a[1]}++;
  $count1++;
}
for my $i(keys %h1)
{
   $h1{$i}/=$count1;
   $h1{$i}=sprintf "%.2f",$h1{$i};
#   print $i." ".$h1{$i}." ";
}
#print "\n";
while(<Tum>)
{
  chomp;
  next if /X/;
  next if /sample/;
  my @a=split;
  my $chr=$a[1];
  $h2{$a[1]}++;
  $count2++;
}
for my $i(keys %h2)
{
   $h2{$i}/=$count2;
   $h2{$i}=sprintf "%.2f",$h2{$i};
 #  print $i." ".$h2{$i}." ";
}
#print "\n";
#print("$count1,$count2\n");

my $mm=1;
my $str1="a=c(";
my $str2="b=c(";
while($mm<=22)
{
   #print "chr.$mm\n";
   if($mm!=22)
   {
   my $now="chr".$mm;
   #print "$now\n";
   if(exists($h1{$now})) {$str1.=$h1{$now}.",";}
   else {$str1.="0,";}
   if(exists($h2{$now})) {$str2.=$h2{$now}.",";}
   else {$str2.="0,";}
   }
   else {
   my $now="chr".$mm;
 #  print "$now\n";
   if(exists($h1{$now})) {$str1.=$h1{$now}."";}
   else {$str1.="0";}
   if(exists($h2{$now})) {$str2.=$h2{$now}."";}
   else {$str2.="0";}
   }
   $mm++;
}
$str1.=")";
$str2.=")";
#print "$str1\n$str2\n";
=cut
chomp(my $out=`pwd`);
$out.="/autoRun.R";
open OUT,">$out" or die $!;
print "$out\n";

print OUT qw(setwd("./"))."\n";
print OUT qw(png(file="Image/wholeChr.png",width=2200,height=1200))."\n";
print OUT qw(file1=read.table("/public/home/chshen/auto_Image_Pro/normal_ratio.txt"))."\n";
print OUT qw(file2=read.table("normal_result"))."\n";
print OUT qw(file3=read.table("tumor_result"))."\n";
print OUT "r1=c(file1[,1])\n";
print OUT "a=c(file2[,3])\n";
print OUT "b=c(file3[,3])\n";
print OUT "m=matrix(c(a,b,r1),3,22,byrow=TRUE)\n";
print OUT qw(x=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))."\n";
print OUT qw(barplot(m,names.arg=x,cex.lab=2,col=c("dodgerblue3","violetred3","forestgreen"),ylab="Ratio",legend=c("normal","tumor","excepted"),args.legend=list(cex=3),ylim=c(0,0.4),cex.axis=2,cex.names=2))."\n";
print OUT "dev.off()\n";

`Rscript $out`;











