#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;

my $refgene = shift;
my $gene = shift;

die "perl $0 <hg19_refGene.txt> <gene>\n" unless ($refgene && $gene);

my %h;
open REF, $refgene or die $!;
while(<REF>){
	chomp;
	my @a = split;
	my $chr = $a[2];
	my $gene_name = $a[12];
	my $gene_spos = $a[4];
	$h{$gene_name} = "$chr\t$gene_spos";
#	print "##$gene_name\t##$chr\t##$gene_spos\n";
	
}
close REF;

open GENE, $gene or die $!;
while(<GENE>){
	chomp;
	next if /GENE_NAME/;
	my @a = split;
	my $gene_name = $a[0];
#	print "$gene_name\n";
	my $sample_num = $a[1];
	if(exists($h{$gene_name}))
	{
	my @gene_info = split /\s+/, $h{$gene_name};
	my $chr = $gene_info[0];
	my ($chr_num) = $chr=~ /chr(\w+)/;
	my $gene_spos = $gene_info[1];
	print "hs$chr_num\t$gene_spos\t$gene_spos\t-$sample_num\n";
	}
	#else {print $gene_name." is not find\n";}
}
close GENE;
