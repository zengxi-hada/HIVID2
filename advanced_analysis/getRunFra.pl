#!usr/bin/perl
use warnings;
use strict;

chomp(my $tumor_all=`less -SN tumor.txt |wc -l`);
chomp(my $normal_all=`less -SN normal.txt |wc -l`);

my @file=qw(Normal_Common.distri Normal_NFRS.distri Normal_Rare.distri Tumor_Common.distri Tumor_NFRS.distri Tumor_Rare.distri);
print("$tumor_all,$normal_all\n");

my @a;
for my $i(@file)
{
   chomp(my $m=`perl -lane 'print "\$F[3]" if(\$F[3]==0)' $i |wc -l`);
   push @a,$m;
}
print "@a\n";
my $c=0;
for my $j(@a)
{
   if($c<=2)
   {$j/=$normal_all;}
   if($c>2) {$j/=$tumor_all;}
   $c++;
   $j=sprintf "%.2f",$j;
}
print "@a\n";
chomp(my $out=`pwd`);
$out.="/autoRunFra.R";
open OUT,">$out" or die $!;
print "$out\n";
print OUT qw{setwd("./Image")}."\n";
print OUT qw(png(file="Fragilr.png",width=800))."\n";
print OUT "x=c($a[0],$a[1],$a[2])\n";
print OUT "y=c($a[3],$a[4],$a[5])\n";
print OUT "z=c(0.2,0.07,0.35)\n";
print OUT "m=matrix(c(x,y,z),3,3,byrow=TRUE)\n";
print OUT qw(l=list("Common","NFRS","RARE"))."\n";
print OUT qw(barplot(width=2,m,ylim=c(0,0.7),ylab="Ratio",names.arg=c("Common","NFRS","RARE"),cex.axis=0.8,cex.names=0.8,args.legend=l,legend.text=c("Normal","Tumor","expected"),beside=T,col=c("blue","purple","orange")))."\n";
print OUT "dev.off()\n";
`Rscript autoRunFra.R`;
#for my $i(@file)
#{
#   print("$i\n");
#}
