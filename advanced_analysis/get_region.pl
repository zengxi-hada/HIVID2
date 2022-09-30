#!usr/bin/perl

use strict;
use warnings;

my $file = shift;
my $file1 = shift;

die "perl $0 <file>\n" unless $file;

open FILE,$file or die $!;
open FILE1,$file1 or die $!;
my %h;
my %h1;
my ($count1,$count2);
while(<FILE>){
	chomp;
  	next if /sample/;
	my @a=split;
	if($a[4]=~/3-UTR/) {$h{'3-UTR'}++;}
        if($a[4]=~/intergenic/) {$h{'intergenic'}++;}
        if($a[4]=~/intro/) {$h{'intronic'}++;}
        if($a[4]=~/5-UTR/) {$h{'5-UTR'}++;}
        if($a[4]=~/CDS/) {$h{'CDS'}++;}
        if($a[4]=~/Promoter/) {$h{'Promoter'}++;}
	$count1++;
}
#print "TUR3:$h{'UTR3'}\ndownstream:$h{'downstream'}\nexonic:$h{'exonic'}\nintergenic:$h{'intergenic'}\nintronic:$h{'intronic'}\nupstream:$h{'upstream'}\n";

while(<FILE1>){
        chomp;
	next if /sample/;
        my @a=split;
        if($a[4]=~/3-UTR/) {$h1{'3-UTR'}++;}
        if($a[4]=~/intergenic/) {$h1{'intergenic'}++;}
        if($a[4]=~/intro/) {$h1{'intronic'}++;}
	if($a[4]=~/5-UTR/) {$h1{'5-UTR'}++;}
    	if($a[4]=~/CDS/) {$h1{'CDS'}++;}
        if($a[4]=~/Promoter/) {$h1{'Promoter'}++;}
 	$count2++;
}
print "$count1,$count2\n";
my @a=qw(3-UTR 5-UTR intergenic intronic CDS Promoter);
for my $i(@a)
{
 #  print "i:$i\n";
   if(exists($h{$i})) {$h{$i}=$h{$i}/$count1;$h{$i}=sprintf "%.2f",$h{$i};}
   else {$h{$i}=0;}
   if(exists($h1{$i})) {$h1{$i}=$h1{$i}/$count2;$h1{$i}=sprintf "%.2f",$h1{$i};}
   else {$h1{$i}=0;}
}
#print "a:@a\n";
for my $m(keys %h1)
{
  print "keys:$m\n";
  if(exists($h1{$m}))
  {print $m.":".$h1{$m}."\n";}
}
for my $m(keys %h)
{
  print "keys:$m\n";
  if(exists($h{$m}))
  {print $m.":".$h{$m}."\n";}
}
#print "TUR3:$h1{'UTR3'}\ndownstream:$h1{'downstream'}\nexonic:$h1{'exonic'}\nintergenic:$h1{'intergenic'}\nintronic:$h1{'intronic'}\nupstream:$h1{'upstream'}\n";
#my (@b1,@b2);


chomp(my $dir=`pwd`);
my $out=$dir."/autoRun.R";
print "$out\n";

#print "$out\n";
open RR,">$out" or die $!;
print RR qw(setwd("./Image"))."\n";
print RR qw(png(file="gene_Region.png"))."\n";
print RR "x=c($h{'3-UTR'},$h1{'3-UTR'})\n";
print RR "y=c($h{'5-UTR'},$h1{'5-UTR'})\n";
print RR "x1=c($h{'intergenic'},$h1{'intergenic'})\n";
print RR "y1=c($h{'intronic'},$h1{'intronic'})\n";
print RR "x2=c($h{'CDS'},$h1{'CDS'})\n";
print RR "y2=c($h{'Promoter'},$h1{'Promoter'})\n";

print RR "m=matrix(c(x,y,x1,y1,x2,y2),2,6)\n";
print RR qw(barplot(width=2,m,ylab="Ratio",ylim=c(0,0.35),names.arg=c("UTR3","UTR5","intergenic","intronic","CDS","Promoter"),cex.axis=0.8,cex.names=0.8,legend.text=c("Normal","Tumor"),beside=T,col=c("blue","purple")))."\n";
print RR "dev.off()\n";


`Rscript $out`;



