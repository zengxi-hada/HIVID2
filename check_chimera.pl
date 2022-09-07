#!/usr/bin/perl

use strict;
use warnings;

my $human_sam = shift;
my $virus_sam = shift;
my $human_initial = shift;
my $virus_initial = shift;
my $use_human_initial = shift;
my $use_virus_initial = shift;
my $human_initial_fusion = shift;
my $virus_initial_fusion = shift;

die "perl $0 <human.sam> <virus.sam> <human.initial> <virus.initial> <out.human.initial> <out.virus.initial> <out.human.fusion.initial> <out.virus.fusion.initial>\n"  unless ($human_sam && $virus_sam && $human_initial && $virus_initial && $use_human_initial && $use_virus_initial && $human_initial_fusion && $virus_initial_fusion);

my %origin_ids;
open HI, $human_initial or die $!;
while(<HI>){
	chomp;
	next if (/reads/);
	my @a = split;
	my $read_id = $a[-1];
	$origin_ids{$read_id} = 1;
}
close HI;

my %h_as;
open HSM, $human_sam or die $!;
while(<HSM>){
	chomp;
	my @a = split;
    if($a[0]=~/@/){next;}
    my $hbv_cigar = $a[5];
    my $mapq = $a[4];
    my $id = $a[0];
    my $hbv_pos = $a[3];
    my $hbv_ref = $a[2];
	my $flag = $a[1];
    my $hbv_fqseq = $a[9];
	if(exists $origin_ids{$id} && ($flag & 2048)){
		next if $h_as{"human"}{$id};
		my ($align_score) = $_=~/AS:i:(\d+)\s+/;
		$h_as{"human"}{$id} = $align_score;
#		print "$id\t$align_score\n";
	}
}
close HSM;

open VSM, $virus_sam or die $!;
while(<VSM>){
    chomp;
    my @a = split;
    if($a[0]=~/@/){next;}
    my $virus_cigar = $a[5];
    my $mapq = $a[4];
    my $id = $a[0];
    my $virus_pos = $a[3];
    my $virus_ref = $a[2];
    my $virus_fqseq = $a[9];
	if(exists $origin_ids{$id}){
		my ($align_score) = $_=~/AS:i:(\d+)\s+/;
		$h_as{"virus"}{$id} = $align_score;
#		print "$id\t$align_score\n";
	}
}
close VSM;

open UHI, ">$use_human_initial" or die $!;
open UVI, ">$use_virus_initial" or die $!;
open FHI, ">$human_initial_fusion" or die $!;
open FVI, ">$virus_initial_fusion" or die $!;


my %record_fusion;
open HI, $human_initial or die $!;
while(<HI>){
    chomp;
    next if (/reads/);
    my @a = split;
    my $read_id = $a[-1];
#	print "#$human_score\t#$virus_score\n";
	if(exists $h_as{"human"}{$read_id}){				# only consider the reads with flag of 2048
		my $human_score =  $h_as{"human"}{$read_id};
	    my $virus_score = $h_as{"virus"}{$read_id};
		if($virus_score - $human_score >= 10){
			print UHI "$_\n";
		}else{
			print FHI "$_\n";
			$record_fusion{$read_id} = 1;				# record the reads which is rearrangement(嵌合 / fusion) in human genome rather than virus integration
		}  
	}else{												# if no flag of 2048 then the read is good and will used later
		print UHI "$_\n";
	}
}
close HI;
close UHI;
close FHI;

open VI, $virus_initial or die $!;
while(<VI>){
	chomp;
	next if (/reads/);
    my @a = split;
    my $read_id = $a[-1];
	if(exists $record_fusion{$read_id}){                # if the reads is rearrangement(嵌合 / fusion) in human genome
		print FVI "$_\n";
    }else{                                              # if no flag of 2048, then the read is good and will used later
        print UVI "$_\n";
    }
}
close VI;
close UVI;
close FVI;


