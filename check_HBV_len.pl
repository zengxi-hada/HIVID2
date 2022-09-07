#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use Cwd qw/abs_path/;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

=head1 function

   this program is designed to check the hbv fragment length and hbv genotype according to the alignment pos of virus sequence on left and right suport    reads

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

my (%hk, $unit);

open HK, $human_bk or die $!;
while(<HK>){
	next if /ref/;
	my @a = split;
	my $ref = $a[0];
 	my $pos = $a[1];
	my $left_support = $a[2];
    my $right_support = $a[3];
    my $sum_support = $a[4];
	my $norm_support = $a[7];
	$unit = $norm_support/$sum_support;
    my $l_reads = $a[-2];
    my $r_reads = $a[-1];
	my (@l_read_array, @r_read_array);
    @l_read_array = split /,/, $l_reads if ($l_reads ne "0");
	@r_read_array = split /,/, $r_reads if ($r_reads ne "0");
#	$hk{$ref}{$pos} = "$left_support\t$right_support\t$sum_support\t$select_read";
	next if @l_read_array == 0;
	next if @r_read_array == 0;
	my (@ltmp, @rtmp);
	for my $i(@l_read_array){
		$i =~ s/repeat_reads_//;
		push @ltmp, $i;
	}
	for my $i(@r_read_array){
        $i =~ s/repeat_reads_//;
        push @rtmp, $i;
    }
	my $l_select_read = join ",", @ltmp;
	my $r_select_read = join ",", @rtmp;
	$hk{$ref}{$pos} = "$l_select_read\t$r_select_read";
}
close HK;

open BK, $hbv_bk or die $!;
open OUT, ">$out" or die $!;

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

print OUT "human_ref\thuman_Pos\thbv_left_pos\thbv_right_pos\thbv_length\n";

for my $ref(keys %hk){
	for my $pos(keys %{$hk{$ref}}){
		my $lr = [split /\t/, $hk{$ref}{$pos}] -> [0];
		my $rr = [split /\t/, $hk{$ref}{$pos}] -> [1];
#		next if not exists $bk{$lr};
#		next if not exists $bk{$rr};
		my ($hbv_lp, $hbv_rp, $hbv_lref, $hbv_rref, $hbv_len, @lr_array, @rr_array, @l_pos, @r_pos);
		@lr_array = split /,/, $lr;
		@rr_array = split /,/, $rr;
		for my $r1(@lr_array){
        	$r1 =~ s/repeat_reads_//;
	        next if not exists $bk{$r1};
   	        for my $g(@{$bk{$r1}}){
        	    push @l_pos, $g;
    	    }
	    }
	    for my $r2(@rr_array){
	        $r2 =~ s/repeat_reads_//;
	        next if not exists $bk{$r2};
   	        for my $g(@{$bk{$r2}}){
	            push @r_pos, $g;
   		    }
	    }
        ## find the HBV pos of the left support reads set ##
		my (%h1);
		for my $j(@l_pos){
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
		## find the HBV pos of the right support reads set ##
		my (%h2);
		for my $j(@r_pos){
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
		##  check the combination of pos between the left and right reads set ##
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
						my @sort = sort {$a<=>$b} ($jp1, $jp2);
    	                my $sort_key = join "\t", @sort;
	                    next if exists $record{"$jf1\t$sort_key"};
	                    next if ($count2 < 2*$unit);
						next if ($jp1 == $jp2);
	                    my $len = abs($jp1 - $jp2);
	                    push @combination_hbv, "$jf1\t$jp1\t$jp2\t$len";
						$record{"$jf1\t$sort_key"} = 1;
	                }
	            }
			}
		}

		## out put ##
		for my $k(@combination_hbv){
     	   print OUT "$ref\t$pos  --  $k\n";
	    }

#		for my $k1(@{$bk{$lr}}){
#			$hbv_lp = [split /\t/, $k1] -> [1];
#			$hbv_lref = [split /\t/, $k1] -> [0];
#			for my $k2(@{$bk{$rr}}){
#				$hbv_rp = [split /\t/, $k2] -> [1];
#				$hbv_rref = [split /\t/, $k2] -> [0];
#				$hbv_len = abs($hbv_rp - $hbv_lp);
#				print OUT "$ref\t$pos --  $hbv_lref\_$hbv_lp\t$hbv_rref\_$hbv_rp\t$hbv_len\n";
#			}
#		}
	}
}

my $count = 0;

close OUT;
