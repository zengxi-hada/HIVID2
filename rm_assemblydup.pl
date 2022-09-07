#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;

my $fq = shift;
my $bk_file = shift;
my $out = shift;
my $record_dup_file = shift;

die "perl $0 <fq file>  <bk_file>  <OUT PUT> <record_dup_reads_file>\n" unless ($fq && $bk_file && $out && $record_dup_file);

if ($fq =~ /gz/){
	open FQ, "<:gzip", $fq or die $!;
}else{
	open FQ, $fq or die $!;
}

my %readsID;

open BK, $bk_file or die $!;
while(<BK>){
	next if /ref/;
	chomp;
	my @a = split;
	my $l_id = $a[-2];
	my $r_id = $a[-1];
	my $id;
	if($r_id eq "0"){
		$id = $l_id;
	}elsif($l_id eq "0"){
		$id = $r_id;
	}else{
		$id = $l_id.",$r_id";
	}
	my @tmp = split /,/, $id;
	for my $i(@tmp){
		$readsID{$i} = 1;
	}
}
close BK;

#$/ = "@";
#<FQ>;
while(<FQ>){
	my $id_str = $_;
    chomp (my $seq = <FQ>);
	<FQ>;
	<FQ>;
#    my @a = split /\n/, $_;
    my $id = (split /\s+/, $id_str)[0];
#    my $seq = $a[1];
	$id =~ s/^\@//;
#	print "$id\t$seq\n";
    if(exists $readsID{$id}){
        $readsID{$id} = $seq;
    }
}
close FQ;

my %record_dup;
my @key = keys %readsID;
my ($seq1, $seq2, $reverse_complementary_seq1);
for my $i1(@key){
	$seq1 = $readsID{$i1};
	$reverse_complementary_seq1 = $seq1;
	$reverse_complementary_seq1 = reverse ($reverse_complementary_seq1);
	$reverse_complementary_seq1 =~ tr/ATCGatcg/TAGCtagc/;
	next if (exists $record_dup{$i1});
	for my $i2(@key){
		next if($i1 eq $i2);
		$seq2 = $readsID{$i2};
		if($seq1 eq $seq2 || $reverse_complementary_seq1 eq $seq2){
#			$record_dup{$i1} = 1;
			$record_dup{$i2} = 1;
#			print "$i1\t$i2\n";
			
		}
	}
#	delete $record_dup{$i};
}


open OUT, ">$out" or die $!;

print OUT "ref\tpos\tleft_support\tright_support\ttotal_support\tnorm_left\tnorm_right\tnorm_sum\tleft_reads_ID\tright_reads_ID\n";

$/ = "\n";
open BK, $bk_file or die $!;
while(<BK>){
#	print "zengxi\n";
    next if /ref/;
    chomp;
    my @a = split;
    my $l_id = $a[-2];
    my $r_id = $a[-1];
 	my (@l_tmp, @r_tmp);
    @l_tmp = split /,/, $l_id if($l_id ne "0");
	@r_tmp = split /,/, $r_id if($r_id ne "0");
	my (@l_print, @r_print, $lp, $rp);
	if($l_id ne "0"){
	    for my $i(@l_tmp){
    	    next if(exists $record_dup{$i});
			push @l_print, $i;
	    }
	}
	if($r_id ne "0"){
        for my $i(@r_tmp){
            next if(exists $record_dup{$i});
            push @r_print, $i;
        }
    }
	my ($left_support,$right_support, $sum_support, $norm_right, $norm_left, $norm_sum);
	$norm_left = "0.000";
	$norm_right = "0.000";
	$norm_sum = "0.000";
	$left_support = @l_print;
	$right_support = @r_print;
	$sum_support = $left_support + $right_support;
	$norm_left = sprintf("%.3f", ($left_support/$a[2])*$a[5]) if ($a[2] != 0);
	$norm_right = sprintf("%.3f", ($right_support/$a[3])*$a[6]) if ($a[3] != 0);
	$norm_sum = sprintf("%.6f", ($sum_support/$a[4])*$a[7]) if ($a[4] != 0);
	$lp = ($l_id eq "0") ? 0 : (join ",", @l_print);
	$rp = ($r_id eq "0") ? 0 : (join ",", @r_print);
	next if $sum_support==0;
	print OUT "$a[0]\t$a[1]\t$left_support\t$right_support\t$sum_support\t$norm_left\t$norm_right\t$norm_sum\t$lp\t$rp\n";	
}
close BK;
close OUT;

open RD, ">$record_dup_file" or die $!;
for my $read_id(keys %record_dup){
	print RD "$read_id\n";
}
close RD;
