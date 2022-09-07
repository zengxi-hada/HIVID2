#!/usr/bin/perl

use strict;
use warnings;

my $infile = shift;
my $stat = shift;
my $outfile = shift;
my $outeff = shift;
my $boole = shift;
my $distance = shift;

die "perl $0 <initialized breakpoint file> <style_box.stat> <out.bk.final> <out file of effected fragment percentage> [\$boole](output effected fragment percentage or not, y or n) <consistency distance>\n" unless ($infile && $outfile && $stat && $outeff && $boole && $distance);

########################################################################################################
#### this program is to merge breakpoints with distance less than 20bps for the virus genome  ##########
###  Author: zeng xi																		  ##########
###	 Last update date: 2020-12																  ##########
########################################################################################################

open IN, $infile or die $!;
open STAT, $stat or die $!;
open OUT, ">$outfile" or die $!;

my %h;
my %idh;
my %support;
my %count_id;
my %left_right;
my %record_pos;
my %repeat;
my @tmp_record;
my $max_item;
my $count_line = 0;
my (%support_final, %idh_final);

while(<IN>){
    next if (/ref/);
	$count_line++;
    my @a = split;
	next if !@a;
    my $ref = $a[0];
    my $pos = $a[1];
    my $left_sup = $a[2];
    my $right_sup = $a[3];
    my $sum_sup = $a[4];
    my $left_normsup = $a[5];
    my $right_normsup = $a[6];
    my $sum_normsup = $a[7];
    my $left_id = $a[-2];
    my $right_id = $a[-1];
    $left_id = ($left_id eq "0") ? "" : $left_id;
    $right_id = ($right_id eq "0") ? "" : $right_id;
    my @l_ids = split /,/, $left_id;
    my @r_ids = split /,/, $right_id;

    my $flag = 0;
    my $flag2 = 0;
    my $flag3 = 0;
    my @ids = (@l_ids, @r_ids);
    for my $tmpi(@ids){$count_id{$tmpi} = 1;}
    for my $i($pos-$distance..$pos+$distance){
        if(exists $record_pos{$ref}{$i}){
            $flag2 = 1;
            push @tmp_record, "$ref\t$pos\t$sum_sup";
			last;
        }
    }
    if($flag2 == 0){
        $record_pos{$ref}{$pos} = 1;
        if(@tmp_record != 0){
            $max_item =  &max_sup (@tmp_record);
            my $iref = (split /\t/, $max_item)[0];
            my $ipos = (split /\t/, $max_item)[1];
            my $isup = (split /\t/, $max_item)[2];
			$support_final{$iref}{$ipos}{"sum"} = $support{$iref}{$ipos}{"sum"};
            $support_final{$iref}{$ipos}{"left"} = $support{$iref}{$ipos}{"left"};
            $support_final{$iref}{$ipos}{"right"} = $support{$iref}{$ipos}{"right"};
            $support_final{$iref}{$ipos}{"normleft"} = $support{$iref}{$ipos}{"normleft"};
            $support_final{$iref}{$ipos}{"normright"} = $support{$iref}{$ipos}{"normright"};
            $support_final{$iref}{$ipos}{"normsum"} = $support{$iref}{$ipos}{"normsum"};
            push @{$idh_final{$iref}{$ipos}{"left"}}, @{$idh{$iref}{$ipos}{"left"}} if (exists $idh{$iref}{$ipos}{"left"});
            push @{$idh_final{$iref}{$ipos}{"right"}}, @{$idh{$iref}{$ipos}{"right"}} if (exists $idh{$iref}{$ipos}{"right"});
            for my $k(@tmp_record){
                my $iref2 = (split /\t/, $k)[0];
                my $ipos2 = (split /\t/, $k)[1];
                my $isup2 = (split /\t/, $k)[2];
#               print "mark2\n" if ($ipos2 =~ /chr/);
                if($k ne $max_item){
                    $support_final{$iref}{$ipos}{"sum"} += $support{$iref2}{$ipos2}{"sum"};
                    $support_final{$iref}{$ipos}{"left"} += $support{$iref2}{$ipos2}{"left"};		# the sup reads will added up within 20 bps
                    $support_final{$iref}{$ipos}{"right"} += $support{$iref2}{$ipos2}{"right"};
                    $support_final{$iref}{$ipos}{"normleft"} += $support{$iref2}{$ipos2}{"normleft"};
                    $support_final{$iref}{$ipos}{"normright"} += $support{$iref2}{$ipos2}{"normright"};
                    $support_final{$iref}{$ipos}{"normsum"} += $support{$iref2}{$ipos2}{"normsum"};
                    push @{$idh_final{$iref}{$ipos}{"left"}}, @{$idh{$iref2}{$ipos2}{"left"}};
                    push @{$idh_final{$iref}{$ipos}{"right"}}, @{$idh{$iref2}{$ipos2}{"right"}};
                }
            }
        }
        @tmp_record = ("$ref\t$pos\t$sum_sup");
    }

    $support{$ref}{$pos}{"left"} = $left_sup;
    $support{$ref}{$pos}{"right"} = $right_sup;
    $support{$ref}{$pos}{"sum"} = $sum_sup;
    $support{$ref}{$pos}{"normleft"} = $left_normsup;
    $support{$ref}{$pos}{"normright"} = $right_normsup;
    $support{$ref}{$pos}{"normsum"} = $sum_normsup;
    push @{$idh{$ref}{$pos}{"left"}}, @l_ids;
    push @{$idh{$ref}{$pos}{"right"}}, @r_ids;
}
close IN;

if($count_line == 0){
	die "\nno virus integration in $infile\n\n";
}

if(!@tmp_record){
$max_item =  &max_sup (@tmp_record);
my $iref = (split /\t/, $max_item)[0];
my $ipos = (split /\t/, $max_item)[1];
my $isup = (split /\t/, $max_item)[2];
$support_final{$iref}{$ipos}{"sum"} = $support{$iref}{$ipos}{"sum"};
$support_final{$iref}{$ipos}{"left"} = $support{$iref}{$ipos}{"left"};
$support_final{$iref}{$ipos}{"right"} = $support{$iref}{$ipos}{"right"};
$support_final{$iref}{$ipos}{"normleft"} = $support{$iref}{$ipos}{"normleft"};
$support_final{$iref}{$ipos}{"normright"} = $support{$iref}{$ipos}{"normright"};
$support_final{$iref}{$ipos}{"normsum"} = $support{$iref}{$ipos}{"normsum"};
#print "mark3\n" if ($ipos =~ /chr/);
push @{$idh_final{$iref}{$ipos}{"left"}}, @{$idh{$iref}{$ipos}{"left"}};
push @{$idh_final{$iref}{$ipos}{"right"}}, @{$idh{$iref}{$ipos}{"right"}};
for my $i(@tmp_record){
    my $iref2 = (split /\t/, $i)[0];
    my $ipos2 = (split /\t/, $i)[1];
    my $isup2 = (split /\t/, $i)[2];
#   print "mark4\n" if ($ipos2 =~ /chr/);
    if($i ne $max_item){
        $support_final{$iref}{$ipos}{"sum"} += $support{$iref2}{$ipos2}{"sum"};
        $support_final{$iref}{$ipos}{"left"} += $support{$iref2}{$ipos2}{"left"};
        $support_final{$iref}{$ipos}{"right"} += $support{$iref2}{$ipos2}{"right"};
        $support_final{$iref}{$ipos}{"normleft"} += $support{$iref2}{$ipos2}{"normleft"};
        $support_final{$iref}{$ipos}{"normright"} += $support{$iref2}{$ipos2}{"normright"};
        $support_final{$iref}{$ipos}{"normsum"} += $support{$iref2}{$ipos2}{"normsum"};
        push @{$idh_final{$iref}{$ipos}{"left"}}, @{$idh{$iref2}{$ipos2}{"left"}};
        push @{$idh_final{$iref}{$ipos}{"right"}}, @{$idh{$iref2}{$ipos2}{"right"}};
    }
}
}


my ($library_name) = $stat=~/step3\/(.*)\//;
my @b=grep {/Total/} <STAT>;
my ($total_reads) = $b[0]=~ /\s+(\d+)\n/;
my $count_of_id = keys %count_id;
my $total_reads2 = $total_reads/(1000*1000);

print OUT "ref\tpos\tleft_support\tright_support\ttotal_support\tnorm_left\tnorm_right\tnorm_sum\tleft_reads_ID\tright_reads_ID\n";

for my $ref(keys %support_final){
    for my $i(keys %{$support_final{$ref}}){
        my $l_reads_id = 0;
        my $r_reads_id = 0;
        $l_reads_id = join ",", @{$idh_final{$ref}{$i}{"left"}} if (exists $idh_final{$ref}{$i}{"left"} && @{$idh_final{$ref}{$i}{"left"}}!=0);
        $r_reads_id = join ",", @{$idh_final{$ref}{$i}{"right"}} if (exists $idh_final{$ref}{$i}{"right"} && @{$idh_final{$ref}{$i}{"right"}}!=0);
        my $print_left = $support_final{$ref}{$i}{"left"} ? $support_final{$ref}{$i}{'left'} : 0;
        my $print_right = $support_final{$ref}{$i}{'right'} ? $support_final{$ref}{$i}{'right'} : 0;
        my $sum = $support_final{$ref}{$i}{"sum"};
        my $norm_print_left = sprintf ("%.3f", $support_final{$ref}{$i}{"normleft"});
		my $norm_print_right = sprintf ("%.3f", $support_final{$ref}{$i}{"normright"});
        my $norm_sum = sprintf ("%.6f", $support_final{$ref}{$i}{"normsum"});
        print OUT "$ref\t$i\t$print_left\t$print_right\t$sum\t$norm_print_left\t$norm_print_right\t$norm_sum\t$l_reads_id\t$r_reads_id\n";
    }
}
close OUT;

##############  subroutine    ##########

sub max_sup{
    my @a = @_;
    my $sup = (split /\t/,$a[0])[2];
    my $ref = (split /\t/,$a[0])[0];
    my $pos = (split /\t/, $a[0])[1];
#   print "sup  $sup\n";
#   print "ref  $ref\n";
#   print "pos  $pos\n";
    for my $i(@a){
        my @b = split /\t/, $i;
        if($b[2]>$sup){
            $sup = $b[2];
            $ref = $b[0];
            $pos = $b[1];
        }
    }
    return "$ref\t$pos\t$sup";
}
