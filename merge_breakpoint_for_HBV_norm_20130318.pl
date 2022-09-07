#!/usr/bin/perl

use strict;
use warnings;

##############################################################################################################################################################
###	This program is to cluster reads coverring the virus integration site and report the breakpoints for the host genome								###### 							
###	Author: zeng xi																																		######
### Last update date: 2020-12																															######
##############################################################################################################################################################

my $infile = shift;
my $stat = shift;
my $outfile = shift;
my $outeff = shift;
my $boole = shift;
my $distance = shift;
my $uniq_rate_file = shift;

die "perl $0 <initialized breakpoint file> <style_box.stat> <out.bk.final> <out file of effected fragment percentage> [\$boole](output effected fragment percentage or not, y or n) <consistency distance> <uniq_rate_file>\n" unless ($infile && $outfile && $stat && $outeff && $boole && $distance && $uniq_rate_file);

open IN, $infile or die $!;
open OUT, ">$outfile" or die $!;
open STAT, $stat or die $!;

my %h;
my %idh;
my %support;
my %count_id;
my %record_orientation;
my %left_right;

if($distance eq "n"){$distance = 0;}   # Take the cigar 1M5S and 1S5M as an example. The integration pos of 1M5S is read_pos_1M5S+1; the integration pos of 1S5M is read_pos_1S5M. Because read_pos_1M5S+1 = read_pos_1S5M, here the consistency distance of different integration sites should be 0, Which means in our results the upstreams(left integration site) and downstreams(right) integration site has the same postion. They all equal to ownstreams(right) integration site postion.

my $uniq_rate;
open URF, $uniq_rate_file or die $!;
while(<URF>){
	chomp;
	my @a = split;
	$uniq_rate = $a[1];
}
close URF;
#print "$uniq_rate\n";

while(<IN>){
	next if (/ref/);
	my @a = split;
	my $ref = $a[0]; 
	my $pos = $a[1];
	my $orientation = $a[2];
	my $cigar = $a[3];
	my $id = $a[-1];
	my $flag = 0;
	$count_id{$id} = 1;
	$record_orientation{$ref}{$pos} = $orientation;
	for my $i($pos-$distance..$pos+$distance){
        if(exists $h{$ref}{$i} && $record_orientation{$ref}{$pos} eq "left_support"){
 #           $h{$ref}{$i}++;
 			($id=~/^left/ || $id=~/^trim_se/) ? ($h{$ref}{$i}+=1) : ($h{$ref}{$i}+=1);
            $flag = 1;
            push @{$idh{$ref}{$i}{"left"}}, $id;
            ($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"left"}+=1) : ($support{$ref}{$i}{"left"}+=1);	## if the paired reads supporting the bk were assembled then supported num should add 1, if noth assembel ,then add 1
        }elsif(exists $h{$ref}{$i} && $record_orientation{$ref}{$pos} eq "right_support"){
#            $h{$ref}{$i}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($h{$ref}{$i}+=1) : ($h{$ref}{$i}+=1);
            $flag = 1;
            push @{$idh{$ref}{$i}{"right"}}, $id;
#            $support{$ref}{$i}{"right"}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"right"}+=1) : ($support{$ref}{$i}{"right"}+=1);
		}elsif((exists $h{$ref}{$i}) && ($record_orientation{$ref}{$pos} eq "left_right") && ($cigar!~/\D(\d+)I/)){
#            $h{$ref}{$i}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($h{$ref}{$i}+=1) : ($h{$ref}{$i}+=1);
            $flag = 1;
            if (not exists $left_right{"$id\t$cigar\t$orientation"}){
                push @{$idh{$ref}{$i}{"right"}}, $id;
#                $support{$ref}{$i}{"right"}++;
				($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"right"}+=1) : ($support{$ref}{$i}{"right"}+=1);
            }else{
                push @{$idh{$ref}{$i}{"left"}}, $id;
#                $support{$ref}{$i}{"left"}++;
				($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"left"}+=1) : ($support{$ref}{$i}{"left"}+=1);
            }
            $left_right{"$id\t$cigar\t$orientation"} = 1;
        }elsif((exists $h{$ref}{$i}) && ($record_orientation{$ref}{$pos} eq "left_right") && ($cigar=~/\D(\d+)I/) && ($1<5)){
#            $h{$ref}{$i}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($h{$ref}{$i}+=1) : ($h{$ref}{$i}+=1);
            $flag = 1;
            if (not exists $left_right{"$id\t$cigar\t$orientation"}){
                push @{$idh{$ref}{$i}{"right"}}, $id;
#                $support{$ref}{$i}{"right"}++;
				($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"right"}+=1) : ($support{$ref}{$i}{"right"}+=1);
            }else{
                push @{$idh{$ref}{$i}{"left"}}, $id;
#                $support{$ref}{$i}{"left"}++;
				($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"left"}+=1) : ($support{$ref}{$i}{"left"}+=1);
            }
            $left_right{"$id\t$cigar\t$orientation"} = 1;
        }elsif((exists $h{$ref}{$i}) && ($record_orientation{$ref}{$pos} eq "left_right") && ($cigar!~/\D(\d+)I/)){
#            $h{$ref}{$i}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($h{$ref}{$i}+=1) : ($h{$ref}{$i}+=1);
            $flag = 1;
            if (not exists $left_right{"$id\t$cigar\t$orientation"}){
                push @{$idh{$ref}{$i}{"right"}}, $id;
#                $support{$ref}{$i}{"right"}++;
				($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"right"}+=1) : ($support{$ref}{$i}{"right"}+=1);
            }else{
                push @{$idh{$ref}{$i}{"left"}}, $id;
#                $support{$ref}{$i}{"left"}++;
				($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"left"}+=1) : ($support{$ref}{$i}{"left"}+=1);
            }
            $left_right{"$id\t$cigar\t$orientation"} = 1;
        }elsif((exists $h{$ref}{$i}) && ($record_orientation{$ref}{$pos} eq "left_right") && ($cigar=~/\D(\d+)I/) && $1>=5){			# there is a HBV insertion inside a read. In this case, this read is both right and left support
#            $h{$ref}{$i}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($h{$ref}{$i}+=1) : ($h{$ref}{$i}+=1);
            $flag = 1;
            push @{$idh{$ref}{$i}{"right"}}, $id;
#            $support{$ref}{$i}{"right"}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"right"}+=1) : ($support{$ref}{$i}{"right"}+=1);
            push @{$idh{$ref}{$i}{"left"}}, $id;
#            $support{$ref}{$i}{"left"}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$i}{"left"}+=1) : ($support{$ref}{$i}{"left"}+=1);
            $left_right{"$id\t$cigar\t$orientation"} = 1;
        }	
    }
	if($flag == 0){								## for the last line
#        $h{$ref}{$pos} = 1;
		($id=~/^left/ || $id=~/^trim_se/) ? ($h{$ref}{$pos}+=1) : ($h{$ref}{$pos}+=1);
#       $record_orientation{$ref}{$pos} = $orientation;
#       push @{$idh{$ref}{$pos}}, $id;
        if($orientation eq "left_support"){
#            $support{$ref}{$pos}{"left"}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$pos}{"left"}+=1) : ($support{$ref}{$pos}{"left"}+=1);
            push @{$idh{$ref}{$pos}{"left"}}, $id;
        }elsif($orientation eq "right_support"){
#            $support{$ref}{$pos}{"right"}++;
			($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$pos}{"right"}+=1) : ($support{$ref}{$pos}{"right"}+=1);
            push @{$idh{$ref}{$pos}{"right"}}, $id;
        }else{
            if(not exists $left_right{"$id\t$cigar\t$orientation"}){
#                $support{$ref}{$pos}{"right"}++;
				($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$pos}{"right"}+=1) : ($support{$ref}{$pos}{"right"}+=1);
                push @{$idh{$ref}{$pos}{"right"}}, $id;
            }else{
#                $support{$ref}{$pos}{"left"}++;
				($id=~/^left/ || $id=~/^trim_se/) ? ($support{$ref}{$pos}{"left"}+=1) : ($support{$ref}{$pos}{"left"}+=1);
                push @{$idh{$ref}{$pos}{"left"}}, $id;
            }
            $left_right{"$id\t$cigar\t$orientation"} = 1;
        }
    }
}
close IN;

#my ($library_name) = $stat=~/step3\/(.*)\//;
#my @b=grep {/Total/} <STAT>;
#my ($total_reads) = $b[0]=~ /\s+(\d+)\n/;
#my $count_of_id = keys %count_id;
#my $total_reads2 = $total_reads/(1000*1000);
#close STAT;

my ($library_name) = $stat=~/step3\/(.*)\//;
open STAT, $stat or die $!;
my @b=grep {/Paired-Total/} <STAT>;
my ($total_pair_reads) = $b[0]=~ /\s+(\d+)\n/;
close STAT;

open STAT, $stat or die $!;
my @c=grep {/Unpaired-Total/} <STAT>;
my ($total_unpair_reads) = $c[0]=~ /\s+(\d+)\n/;
close STAT;

my $count_of_id = keys %count_id;
my $total_reads = $total_pair_reads + $total_unpair_reads;
my $total_reads2 = ($total_pair_reads + $total_unpair_reads) * $uniq_rate / (1000*1000);

print OUT "ref\tpos\tleft_support\tright_support\ttotal_support\tnorm_left\tnorm_right\tnorm_sum\tleft_reads_ID\tright_reads_ID\n";

for my $ref(keys %h){
    for my $i(keys %{$h{$ref}}){
        my $l_reads_id = 0;
        my $r_reads_id = 0;
        $l_reads_id = join ",", @{$idh{$ref}{$i}{"left"}} if exists $idh{$ref}{$i}{"left"};
        $r_reads_id = join ",", @{$idh{$ref}{$i}{"right"}} if exists $idh{$ref}{$i}{"right"};
        my $print_left = $support{$ref}{$i}{"left"} ? $support{$ref}{$i}{'left'} : 0;			# 实际的左支持reads数目
        my $print_right = $support{$ref}{$i}{'right'} ? $support{$ref}{$i}{'right'} : 0;		# 实际的右支持reads数目
#        my $sum = $h{$ref}{$i};
		my $sum = $h{$ref}{$i};
        my $norm_print_left = sprintf ("%.3f", $print_left/$total_reads2);
        my $norm_print_right = sprintf ("%.3f", $print_right/$total_reads2);
        my $norm_sum = sprintf ("%.6f", $sum/$total_reads2);
        print OUT "$ref\t$i\t$print_left\t$print_right\t$sum\t$norm_print_left\t$norm_print_right\t$norm_sum\t$l_reads_id\t$r_reads_id\n";
    }
}
close OUT;

if($boole eq "y"){

    my $count_of_id = keys %count_id;
#   my $effect_fragment_percent = sprintf ("%.8f", $count_of_id/$total_reads);
    my $effect_fragment_percent = $count_of_id/$total_reads;

    open OUTEFF, ">$outeff" or die $!;
    print OUTEFF "$library_name\t$effect_fragment_percent\n";
}
