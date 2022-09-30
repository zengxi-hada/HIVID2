#!/usr/bin/perl

use strict;
use warnings;

my $anno = shift;
my $chr_len = shift;

die "perl $0 <anno> <chr_len>\n" unless ($anno && $chr_len);
=head
my %ratio = ('chr1' => 0.00149,
       'chr10' => 0.00149,
        'chr11' => 0.00149,
        'chr12' => 0.00149,
        'chr13' => 0.00149,
        'chr14' => 0.00149,
        'chr15' => 0.00149,
        'chr16' => 0.00149,
        'chr17' => 0.00149,
        'chr18' => 0.00149,
        'chr19' => 0.00149,
        'chr2' => 0.00149,
        'chr20' => 0.00149,
        'chr21' => 0.00149,
        'chr22' => 0.00149,
        'chr3' => 0.00149,
        'chr4' => 0.00149,
        'chr5' => 0.00149,
        'chr6' => 0.00149,
        'chr7' => 0.00149,
        'chr8' => 0.00149,
        'chr9' => 0.00149,
);
=cut
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
my @chr_array = (1..22);

my %len;
open CL, $chr_len or die $!;
while(<CL>){
	next if /Total/;
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $length = $a[1];
	my $spos1 = $length - 2000000;
	my $epos1 = $length;
	my $spos2 = 0;
	my $epos2 = 2000000;
	push @{$len{$chr}}, "$spos1\t$epos1";
	push @{$len{$chr}}, "$spos2\t$epos2";
}
close CL;

my %h;
open ANNO, $anno or die $!;
while(<ANNO>){
	chomp;
        next if /sample/;
	next if /chrM/;
	next if /chrX/;
	next if /chrY/;
	next if /chr23/;
	next if /chr24/;
	my @a = split;
	my $lib = $a[0];
	my $chr = $a[1];
	my $pos = $a[2];
	$h{$lib}{'count'}++;
	for my $pos_str(@{$len{$chr}}){
		my $spos = (split /\t/, $pos_str)[0];
		my $epos = (split /\t/, $pos_str)[1];
		if($pos>=$spos && $pos<=$epos){
#			$h{$lib}{'count'}++;
			$h{$lib}{$chr}++;
		}
	}
}
close ANNO;

for my $i(@chr_array, 'count'){
	for my $lib(keys %h){
		if($i eq "count"){
			if(not exists $h{$lib}{$i}){
				$h{$lib}{$i} = 0;
			}
		}else{
			if(not exists $h{$lib}{"chr$i"}){
				$h{$lib}{"chr$i"} = 0;
			}
		}
	}
}

for my $lib(keys %h){
	print "\t";
	for my $chr(sort keys %{$h{$lib}}){
		print "$chr\t$chr\t"
	}
	print "\n";
	last;
}

#print "lib\texonic\te_exonic\tintergenic\te_intergenic\tintronic\te_intronic\tpromoter\te_promoter\n";
for my $lib(keys %h){
	print "$lib\t";
	my $lib_bknum = $h{$lib}{'count'};
	for my $chr(sort keys %{$h{$lib}}){
		next if($chr eq 'count');
		my $num = $h{$lib}{$chr};
		my $expect =$ratio{$chr};
#		print "$lib_bknum#yangjie\n";
#		print "$ratio{$chr}#zengxi#$chr\n";
		print "$num\t$expect\t";
	}
	print "\n";
}
