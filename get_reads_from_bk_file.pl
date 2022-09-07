#!/usr/bin/perl

use strict;
use warnings;

my $fq = shift;
my $bk_file = shift;
my $target_bk = shift;			# chr:pos

die "pelr $0 <fq.file> <bk_file> <target_bk>\n" unless ($fq && $bk_file && $target_bk);

#my @array_target = split /:/, $target_bk;
#my $target_bk_chr = $array_target[0];
#my $target_bk_pos = $array_target[1];

my %bk_reads;
open BK, $bk_file or die $!;
while(<BK>){
	chomp;
	my @a = split;
	my $chr = $a[0];
    my $pos = $a[1];
    my $lsup = $a[2];
    my $rsup = $a[3];
    my $total_sup = $a[4];
    my $lnorm = $a[5];
    my $rnorm = $a[6];
    my $total_norm = $a[7];
    my $l_reads = $a[-2];
    my $r_reads = $a[-1];
	my @array_l_reads = ($l_reads ne 0) ? (split /,/, $l_reads) : ();
	my @array_r_reads = ($r_reads ne 0) ? (split /,/, $r_reads) : ();
	@{$bk_reads{"$chr:$pos"}} = (@array_l_reads, @array_r_reads);
}

open FQ, $fq or die $!;
while(<FQ>){
	chomp (my $id = $_);
	$id =~ s/\@//;
	chomp (my $seq = <FQ>);
	my $plus = <FQ>;
	my $qual = <FQ>;
#	print "$id\n";
	if(exists $bk_reads{$target_bk}){
		my @reads_id_arr = @{$bk_reads{$target_bk}};
		for my $read_id(@reads_id_arr){
			if($read_id eq $id){
				print "$read_id\t$seq\n";
			}
		}
	}
}
close FQ;





