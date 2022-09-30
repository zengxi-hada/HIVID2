#!usr/bin/perl
use warnings;
use strict;

chomp(my $dir=`pwd`);
#print "$dir\n";
#my $tumor=$dir."/tumor.bk.cpg.distr |wc -l";
chomp(my $a=`perl -lane 'print "\$F[3]" if(\$F[3]==0)' $dir/tumor.bk.cpg.distr |wc -l`);
my $normal=$dir."/normal.bk.cpg.distr |wc -l";
chomp(my $b=`perl -lane 'print "\$F[3]" if(\$F[3]==0)' $normal`);
#print "$normal,$tumor\n";
print "$a,$b\n";
chomp(my $normal_all=`less -SN $dir/normal.bk.cpg.distr |wc -l`);
chomp(my $tumor_all=`less -SN $dir/tumor.bk.cpg.distr|wc -l`);
#eval{
	$a/=$normal_all;
	$b/=$tumor_all;
#}
$a=sprintf "%.2f",$a;
$b=sprintf "%.2f",$b;
my $out=$dir."/RunCPG.R";
print "$out\n";
open RR,">$out" or die $!;
print RR qw(setwd("./"))."\n";
print RR "png(file=\"Image/CPG.png\")\n";
print RR "x=c($a,$b,0.1)\n";
print RR "barplot(x,ylab=\"ratio\",ylim=c(0,0.7),names.arg=c(\"normal\",\"tumor\",\"except\"),col=c(\"blue\",\"purple\",\"yellow\"),main=\"CPG\")\n";
print RR "dev.off()\n";
#setwd("./")
#png(file="result1.png")
#x=c(0.43,0.57,0.1)
#barplot(x,ylab="ratio",ylim=c(0,0.7),names.arg=c("normal","tumor","except"),col=c("blue","purple","yellow"),main="CPG")
#dev.off()
`Rscript RunCPG.R`;
