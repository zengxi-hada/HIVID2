#!/usr/bin/perl

use strict;
use warnings;

################################################################################################################################################################################
#### This program is the first round of PCR duplication removement for the supporting reads of breakpoints																########
#### Author: zeng xi																																					########
#### last Update: 2020-12																																				########
################################################################################################################################################################################################################################################################################################################################################################
#### This program is the first round of PCR duplication removement for the supporting reads of breakpoints																########
#### Author: zeng xi																																					########
#### last Update: 2020-12																																				########
################################################################################################################################################################################################################################################################################################################################################################
#### This program is the first round of PCR duplication removement for the supporting reads of breakpoints																########
#### Author: zeng xi																																					########
#### last Update: 2020-12																																				########
################################################################################################################################################################################

my $human_bk = shift;
my $virus_bk = shift;
my $human_sam = shift;
my $virus_sam = shift;
my $human_bk_rmdup = shift;
my $virus_bk_rmdup = shift;

die "perl $0 <human.bk/tumor.bk> <virus.bk/normal.bk> <human_sam> <virus_sam> <human_bk.rmdup/tumor.bk.rmdup> <virus_bk.rmdup/normal_bk.rmdup>\n" unless ($human_bk && $virus_bk && $virus_sam && $human_sam && $human_bk_rmdup && $virus_bk_rmdup);

my %h_hbs;
my %h_hbs_reverse;
open HBS, $virus_sam or die $!;
while(<HBS>){
	next if /\@SQ/;
	next if /\@PG/;
	chomp;
	my @a = split;
#	my $lib = $a[0];
	next if $a[5]=~/H/;
	my $read_id = $a[0];
	my $chr = $a[2];
    my $pos = $a[3];
    my $cigar = $a[5];
	my $len = length($a[9]);
	$h_hbs{$read_id} = "$chr\t$pos\t$len";
#	print "$read_id\n";
	push @{$h_hbs_reverse{"$chr\t$pos\t$len"}}, $read_id;
}
close HBS;

my %h_hus;
my %h_hus_reverse;
open HUS, $human_sam or die $!;
while(<HUS>){
	next if /\@SQ/;
	next if /\@PG/;
	chomp;
	my @a = split;
	next if $a[5]=~/H/;				# if there is hard clip, pass the reads
#	my $lib = $a[0];
	my $read_id = $a[0];
	my $chr = $a[2];
	my $pos = $a[3];
	my $cigar = $a[5];
	my $len = length($a[9]);
	$h_hus{$read_id} = "$chr\t$pos\t$len";
#	print "$read_id\n";
	push @{$h_hus_reverse{"$chr\t$pos\t$len"}}, $read_id;
}
close HUS;

my %record_id;
my %record_uniq_id;
open HUBK, $human_bk or die $!;
open OHUBK, ">$human_bk_rmdup" or die $!;
print OHUBK "ref\tpos\tleft_support\tright_support\ttotal_support\tnorm_left\tnorm_right\tnorm_sum\tleft_reads_ID\tright_reads_ID\n";

while(<HUBK>){
	next if /left_reads_ID/;
	chomp;
	my @a = split;
#	my $lib = $a[1];
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
	my ($l_uniq_read_array_id_str, $r_uniq_read_array_id_str, $l_count_rmdup, $r_count_rmdup);
	if($l_reads eq "0"){
		$l_uniq_read_array_id_str = 0;
		$l_count_rmdup = 0;
	}else{		
	    my @l_read_array_id = split /,/, $l_reads;
		my @l_read_array_cigar = map {$h_hus{$_}} @l_read_array_id;
		for my $i(@l_read_array_id){
			$record_id{$h_hus{$i}}=$i;
		}
		my %count;
#		print "@l_read_array_id\n";
#		print "@l_read_array_cigar\n";
		my @l_uniq_read_array_cigar = grep { ++$count{$_} < 2 } @l_read_array_cigar;		# if the alignment position and aligned len are the same,then the two reads is considered as PCR dup
###		print "@l_uniq_read_array_cigar\n";
#		my @l_uniq_read_array_id = map {my $human_cigar=(split /,/, $_)[1]; ${$h_hus_reverse{$human_cigar}}[0]} @l_uniq_read_array_cigar;
		my @l_uniq_read_array_id = map {my $read_id=$record_id{$_}; $read_id} @l_uniq_read_array_cigar;

		for my $uniq_id(@l_uniq_read_array_id){ $record_uniq_id{$uniq_id} = 1; }

#		print "zengxi\t@l_uniq_read_array_cigar\n";
		$l_uniq_read_array_id_str = join (",", @l_uniq_read_array_id);
#		$l_count_rmdup = @l_uniq_read_array_id;
		for my $xx(@l_uniq_read_array_id){
            if($xx=~/^left/ || $xx=~/^trim_se/){
                $l_count_rmdup += 1;
            }else{
                $l_count_rmdup += 1;
            }
        }
	}

	if($r_reads eq "0"){
        $r_uniq_read_array_id_str = 0;
		$r_count_rmdup = 0;
    }else{
        my @r_read_array_id = split /,/, $r_reads;
        my @r_read_array_cigar = map {$h_hus{$_}} @r_read_array_id;
		for my $i(@r_read_array_id){
            $record_id{$h_hus{$i}}=$i;							# store the pos and cigar and read id 
        }
        my %count;
        my @r_uniq_read_array_cigar = grep { ++$count{$_} < 2 } @r_read_array_cigar;		# if the alignment position and aligned len are the same,then the two reads is considered as PCR dup
##        my @r_uniq_read_array_id = map {my $human_cigar=(split /,/, $_)[1]; ${$h_hus_reverse{$human_cigar}}[0]} @r_uniq_read_array_cigar;
		 my @r_uniq_read_array_id = map {my $read_id=$record_id{$_}; $read_id} @r_uniq_read_array_cigar;
 
		for my $uniq_id(@r_uniq_read_array_id){ $record_uniq_id{$uniq_id} = 1; }

        $r_uniq_read_array_id_str = join (",", @r_uniq_read_array_id);
#		$r_count_rmdup = @r_uniq_read_array_id;
		for my $xx(@r_uniq_read_array_id){
			if($xx=~/^left/ || $xx=~/^trim_se/){
				$r_count_rmdup += 1;
			}else{
				$r_count_rmdup += 1;
			}
		}	
    }
	pop @a; pop @a;
	my $print_str = join ("\t", @a);
#	print OHUBK "$print_str\t\t$l_count_rmdup\t$r_count_rmdup\t$l_uniq_read_array_id_str\t$r_uniq_read_array_id_str\n";
	$l_reads = ($l_reads eq "0") ? "" : $l_reads;
	my $total_count_rmdup = $l_count_rmdup + $r_count_rmdup;
	my $new_lnorm = ($lsup eq "0") ? "0.000" : sprintf("%.3f", $lnorm*$l_count_rmdup/$lsup);
	my $new_rnorm = ($rsup eq "0") ? "0.000" : sprintf("%.3f", $rnorm*$r_count_rmdup/$rsup);
	my $new_total_norm = sprintf("%.3f", $total_norm * $total_count_rmdup/$total_sup);
#	print "$lsup\t$lnorm\t$new_lnorm\n";
#	print "$total_sup\t$total_count_rmdup\t$new_total_norm\n";
    print OHUBK "$chr\t$pos\t$l_count_rmdup\t$r_count_rmdup\t$total_count_rmdup\t$new_lnorm\t$new_rnorm\t$new_total_norm\t$l_uniq_read_array_id_str\t$r_uniq_read_array_id_str\n";
}
close HUBK;
close OHUBK;

open VBK, $virus_bk or die $!;
open OVBK, ">$virus_bk_rmdup" or die $!;
print OVBK "ref\tpos\tleft_support\tright_support\ttotal_support\tnorm_left\tnorm_right\tnorm_sum\tleft_reads_ID\tright_reads_ID\n";
while(<VBK>){
	next if /left_reads_ID/;
    chomp;
    my @a = split;
#    my $lib = $a[1];
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
	
	my @l_tmp = split /,/, $l_reads if($l_reads ne "0");
    my @r_tmp = split /,/, $r_reads if($r_reads ne "0");
    my (@l_print, @r_print, $lp, $rp, $left_support, $right_support);
    if($l_reads ne "0"){
        for my $i(@l_tmp){
            if(exists $record_uniq_id{$i}){
            	push @l_print, $i;
				if($i=~/^left/ || $i=~/^trim_se/){
					$left_support += 1;
				}else{
					$left_support += 1;
				}
			}
        }
    }
    if($r_reads ne "0"){
        for my $i(@r_tmp){
            if(exists $record_uniq_id{$i}){
            	push @r_print, $i;
				if($i=~/^left/ || $i=~/^trim_se/){
                    $right_support += 1;
                }else{
                    $right_support += 1;
                }
			}
        }
    }
    my ($sum_support, $norm_right, $norm_left, $norm_sum);
    $norm_left = "0.000";
    $norm_right = "0.000";
    $norm_sum = "0.000";
    $left_support = 0 if(($l_reads eq "0") || (@l_print==0));
    $right_support = 0 if(($r_reads eq "0") || (@r_print==0));
    $sum_support = $left_support + $right_support;
    $norm_left = sprintf("%.3f", ($left_support/$lsup)*$lnorm) if ($lsup != 0);
    $norm_right = sprintf("%.3f", ($right_support/$rsup)*$rnorm) if ($rsup != 0);
    $norm_sum = sprintf("%.3f", ($sum_support/$total_sup)*$total_norm) if ($total_sup != 0);
    $lp = ($l_reads eq "0") ? 0 : (join ",", @l_print);
    $rp = ($r_reads eq "0") ? 0 : (join ",", @r_print);
    next if $sum_support==0;
    print OVBK "$chr\t$pos\t$left_support\t$right_support\t$sum_support\t$norm_left\t$norm_right\t$norm_sum\t$lp\t$rp\n";

#    print OVBK "$chr\t$pos\t$l_count_rmdup\t$r_count_rmdup\t$total_count\t$new_lnorm\t$new_rnorm\t$new_total_norm\t$l_uniq_read_array_id_str\t$r_uniq_read_array_id_str\n";
}
close VBK;
close OVBK;

################################# subroutine #################################
sub get_uniq{
	my @a = @_;
	for my $i(@a){
		my @b = split /\t/, $i;
		my $lib = $b[0];
		my $chr = $b[1];
		my $pos = $b[2];
		my $cigar = $b[3];
	}
}
