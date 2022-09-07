#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use Cwd qw/abs_path/;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;


=head1 function

   this program is designed to check if there exists replacement between HBV genome and human genome

=head1 usage
      perl $0
			-hk	<str>		the path of the human breakpoint file
			-bk	<str>		the path of the HBV breakpoint file
			-o	<str>		the output path

=head1 version
   v1.0

=head1 author
   zengxi: zengxi@genomics.cn

=head1 date
   2012-02-19

=cut

my ($human_bk, $hbv_bk, $out);

GetOptions(
	'hk=s' => \$human_bk,
	'bk=s' => \$hbv_bk,
	'o=s'  => \$out,
);
	
die `pod2text $0` unless ($human_bk && $hbv_bk && $out);

my %hk;

my ($cutoff, $unit);

open HK, $human_bk or die $!;
while(<HK>){
	next if /ref/;
	my @a = split;
	my $ref = $a[0];
 	my $pos = $a[1];
	my $left_support = $a[2];
    my $right_support = $a[3];
    my $sum_support = $a[4];
	my $norm_sup = $a[7];
	$unit = $norm_sup/$sum_support;
    my $l_reads = $a[-2];
	my $r_reads = $a[-1];
	$l_reads = ($l_reads eq "0") ? "" : $l_reads;
	$r_reads = ($r_reads eq "0") ? "" : $r_reads; 
	my $reads = ($l_reads eq "") ? $r_reads : $l_reads.",$r_reads";
    my @read_array = split /,/, $reads;
#    my $select_read = join ",", @read_array;
	my @tmp_array;
	for my $i(@read_array){
		if($i=~/repeat_reads_/){$i =~ s/repeat_reads_//;}
		push @tmp_array, $i;
	}
#	print Dumper(@read_array),"\n";
	my $select_read = join ",", @tmp_array;
	$hk{$ref}{$pos} = "$left_support\t$right_support\t$sum_support\t$norm_sup\t$select_read";
}
close HK;

open HK, $human_bk or die $!;
open BK, $hbv_bk or die $!;
open OUT, ">$out" or die $!;

my %replacement;
my %re_dup;

while(<HK>){
	next if (/ref/);
	my @a = split;
	my $ref = $a[0];
	my $pos = $a[1];
	my $left_support = $a[2];
	my $right_support = $a[3];
	my $sum_support = $a[4];
	my $norm_sup = $a[7];
	my $l_reads = $a[-2];
    my $r_reads = $a[-1];
	$l_reads = ($l_reads eq "0") ? "" : $l_reads;
    $r_reads = ($r_reads eq "0") ? "" : $r_reads;
	my $reads = ($l_reads eq "") ? $r_reads : $l_reads.",$r_reads";
	my @read_array = split /,/, $reads;
#	my $select_read = $read_array[0];
	my @tmp_array;
	for my $i(@read_array){
        if($i=~/repeat_reads_/){$i =~ s/repeat_reads_//;}
        push @tmp_array, $i;
    }
	my $select_read = join ",", @tmp_array;
	my $flag = 0;
	if(exists $hk{$ref}){
		for my $i(keys %{$hk{$ref}}){
			next if $i == $pos;
			my $tmp = $hk{$ref}{$i};
			my @b = split /\t/, $tmp;
			my $s_read = $b[-1];
			my $lsup = $b[0];
			my $rsup = $b[1];
			my $ssup = $b[2];
			my $nsup = $b[3];
			my @c = sort {$a<=>$b} ($i, $pos);
            my $c = join "\t", @c;
			next if exists $re_dup{"$ref\t$c"};
			next if $select_read eq $s_read;
			if(abs($i - $pos) <= 500){
				my $human_len = abs($i - $pos);
				$replacement{"$select_read\t$s_read"} = "$ref\t$pos\t$left_support\t$right_support\t$sum_support\t$norm_sup  :  $i\t$lsup\t$rsup\t$ssup\t$nsup\t$human_len";
#				$replacement{"$select_read\t$s_read"} = "$ref\t$pos\t$left_support\t$right_support  :  $i\t$lsup\t$rsup";	
			}
			$re_dup{"$ref\t$c"} = 1;
		}
	}
			
}


close HK;

my %bk;

while(<BK>){
	next if /ref/;
    my @a = split;
    my $ref = $a[0];
    my $pos = $a[1];
    my $sum_support = $a[2];
    my $reads = $a[-1];
    my @read_array = split /,/, $reads;
	for my $i(@read_array){
		$i =~ s/repeat_reads_//;
		push @{$bk{$i}}, "$ref\t$pos\t$sum_support";
	}
}

print OUT "human ref\thuman_pos1\thuman_left_sup1\thuman_right_sup1\thuman_sum_sup1\tnorm_sup  :  human_pos2\thuman_left_sup2\thuman_right_sup2\thuman_sum_sup2\tnorm_sup\thuman_length  --  HBV_ref\tHBV_pos\tHBV_sum_sup\tHBV_length\n";
#print OUT "human ref\thuman_pos1\thuman_left_sup1\thuman_right_sup1  :  human_pos2\thuman_left_sup2\thuman_right_sup2  --  HBV_ref\tHBV_pos\tHBV_sum_sup\n";
for my $i(keys %replacement){
	my @read = split /\t/, $i;
#	print "@read\n";
	my ($p1, $p2, $ref1, $ref2, $hash_value1, $hash_value2, @read1, @read2, @r1_pos, @r2_pos);
	@read1 = split /,/, $read[0];
	@read2 = split /,/, $read[1];
	for my $r1(@read1){
		$r1 =~ s/repeat_reads_//;
		next if not exists $bk{$r1};
		for my $g(@{$bk{$r1}}){
			push @r1_pos, $g;
		}
	}
	for my $r2(@read2){
		$r2 =~ s/repeat_reads_//;
		next if not exists $bk{$r2};
		for my $g(@{$bk{$r2}}){
			push @r2_pos, $g;
		}
	}
	## find the HBV pos of the first reads set ##
	my (%h1);
	for my $j(@r1_pos){
		my $flag = 0;
		my $ipos = [split /\t/, $j] -> [1];
		my $iref = [split /\t/, $j] -> [0];
	    for my $i($ipos-0..$ipos+0){
    	    if(exists $h1{$iref}{$i}){
	            $h1{$iref}{$i}++;
	            $flag = 1;
	        }
	    }
	    if($flag == 0){
	        $h1{$iref}{$ipos} = 1;
	    }
	}
	## find the HBV pos of the second reads set ##
	my (%h2);
    for my $j(@r2_pos){
        my $flag = 0;
        my $ipos = [split /\t/, $j] -> [1];
        my $iref = [split /\t/, $j] -> [0];
        for my $i($ipos-0..$ipos+0){
            if(exists $h2{$iref}{$i}){
                $h2{$iref}{$i}++;
                $flag = 1;
            }
        }
        if($flag == 0){
            $h2{$iref}{$ipos} = 1;
        }
    }
	## 	check the combination of pos between the two reads set ##
	my (@combination_hbv, %record);	
	for my $jf1(keys %h1){
		for my $jp1(keys %{$h1{$jf1}}){
			my $count1 = $h1{$jf1}{$jp1};
			$count1 = $unit*$count1;
			next if ($count1 < 2*$unit);
			for my $jf2(keys %h2){
				next if ($jf1 ne $jf2);
				for my $jp2(keys %{$h2{$jf2}}){
					my $count2 = $h2{$jf2}{$jp2};
					$count2 = $unit*$count2;
					next if ($count2 < 2*$unit);
					my $len = abs($jp1 - $jp2);
					my @sort = sort {$a<=>$b} ($jp1, $jp2);
					my $sort_key = join "\t", @sort;
					next if exists $record{"$jf1\t$sort_key"};
#					push @combination_hbv, "$jf1\t$jp1\t$jp2\t$len" if (abs($jp1-$jp2) < 500 && $jp1!=$jp2);
					push @combination_hbv, "$jf1\t$jp1\t$jp2\t$len" if ($jp1!=$jp2);
					$record{"$jf1\t$sort_key"} = 1;
				}
			}
		}
	}
	## out put ##
	for my $k(@combination_hbv){
		print OUT "$replacement{$i}  --  $k\n";
	}
#	for my $k1(@{$bk{$read[0]}}){
#		$p1 = [split /\t/, $k1] -> [1];
#		$ref1 = [split /\t/, $k1] -> [0];
#		for my $k2(@{$bk{$read[1]}}){
#			$p2 = [split /\t/, $k2] -> [1];
#			$ref2 = [split /\t/, $k2] -> [0];
#			next if $ref1 ne $ref2;
#			if(abs($p1-$p2) < 500 && $p1 != $p2){
#				my $hbv_len = abs($p1 - $p2);
#				print OUT "$replacement{$i}  --  $k1\t$k2\t$hbv_len\n";	
#			}
#		}
#	}
}

close OUT;
