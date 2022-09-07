#!/usr/bin/perl

use strict;
use warnings;

##################################################################################################
##### this program is to find the reads which can both mapped onto human and virus genome #########
##### Author: zeng xi                                                                     #########
##### Last update: 2020-12                                                                #########
###################################################################################################

my $hbv_reads = shift;
my $human_reads = shift;
my $hbv_out = shift;
my $human_out = shift;

die "perl $0 <hbv sam file> <human sam file> <hbv_out> <human_out> \n" unless ($hbv_reads && $human_reads && $hbv_out && $human_out);

open HBV, $hbv_reads or die $!;

my %mutual;

my %count1;
my %tmp_hash;
#my @HBV_type = (a..h); 
my %hbv_type;

while(<HBV>){
	chomp;
	my @a = split;
	if($a[0]=~/@/){next;}
	my $hbv_cigar = $a[5];
	my $mapq = $a[4];
	my $id = $a[0];
	my $hbv_pos = $a[3];
	my $hbv_ref = $a[2]; 
	if($hbv_cigar =~ /\*/){next;}
#	my ($type) = $hbv_ref =~ /\w+\_(\w)\d$/;
	my $type = $hbv_ref;                                     ## test at 10:21 2012-3-27
	my $info = "$hbv_cigar\t$mapq\t$hbv_ref\t$hbv_pos";
	if(not exists $mutual{$id}){
		$mutual{$id} = $info;
#		$hbv_type{$id}{$type} = 1;
		$count1{$id} = 1;
	}else{
		if($mapq > (split /\s+/,$mutual{$id})[1]){
			$mutual{$id} = $info;
#			$hbv_type{$id}{$type} = 1;
		}elsif($mapq == (split /\s+/,$mutual{$id})[1]){
#			if(not exists $hbv_type{$id}{$type}){
				$mutual{"$id\_$count1{$id}"} = $info;
				push @{$tmp_hash{$id}}, "$id\_$count1{$id}";
				$count1{$id}++;
#				$hbv_type{$id}{$type} = 1;
#			}
		}
	}
}

close HBV;


open HUMAN, $human_reads or die $!;

#print "readsID\thuman_ref\thuman cigar\thuman position\tHBV_ref\tHBV cigar\tHBV position\n";

my %print_hbv;
my %print_human;

my %count2;       ## modify at 14:52 2011-11-06
my %record_human;

while(<HUMAN>){
	chomp;
	my @a = split;
	if($a[0]=~/@/){next;}
	my $id = $a[0];
	if($a[5] =~ /\*/){next;}
	if(exists $mutual{$id}){
		my $human_cigar = $a[5];
		my $hbv_cigar = (split /\s+/, $mutual{$id})[0];
		my $hbv_pos = (split /\s+/, $mutual{$id})[3];
		my $hbv_ref = (split /\s+/, $mutual{$id})[2];
		my $human_mapq = $a[4];
		my $human_pos = $a[3];
		my $human_ref = $a[2];
		if(exists $tmp_hash{$id}){
			$print_hbv{$id} = "repeat_reads_$id\t$hbv_ref\t$hbv_cigar\t$hbv_pos";
			for my $i(@{$tmp_hash{$id}}){
				my $hbv_ref_1 = (split /\s+/, $mutual{$i})[2];                                  ## modify at 14:52 2011-11-06
				my $hbv_cigar_1 = (split /\s+/, $mutual{$i})[0];								## modify at 14:52 2011-11-06
				my $hbv_pos_1 = (split /\s+/, $mutual{$i})[3];									## modify at 14:52 2011-11-06
				$print_hbv{$i} = "repeat_reads_$id\t$hbv_ref_1\t$hbv_cigar_1\t$hbv_pos_1";					## modify at 14:52 2011-11-06; at 11:53 2012-01-10
			}
		}else{
			$print_hbv{$id} = "$id\t$hbv_ref\t$hbv_cigar\t$hbv_pos";
		}

		if(not exists $print_human{$id}){
#			$print_human{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$hbv_ref\t$hbv_cigar\t$hbv_pos\t$human_mapq";
			$print_human{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";
			$count2{$id} = 1;
		}else{
			if($human_mapq > (split /\s+/, $print_human{$id})[4]){
				$print_human{$id} = "$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";
			}elsif($human_mapq == (split /\s+/, $print_human{$id})[4]){
				if($count2{$id} == 1){
#					my $id_tmp = (split /\t/, $print_human{$id})[0];
#					my $id_printtmp = "repeat_reads_$id_tmp";
					my @tmp_item = (split /\t/, $print_human{$id}); 	
					$print_human{"$id\_0"} = "repeat_reads_$id\t$tmp_item[1]\t$tmp_item[2]\t$tmp_item[3]\t$tmp_item[4]";
#					my $count_tmp_item = @tmp_item;
#					print $count_tmp_item,"\n";
				}
				$print_human{"$id\_$count2{$id}"} = "repeat_reads_$id\t$human_ref\t$human_cigar\t$human_pos\t$human_mapq";		## modify at 14:52 2011-11-06; at 11:55 2012-01-10
				$count2{$id}++;
				$record_human{$id} = 1;
			}
		}
	}
}

close HUMAN;

for my $i(keys %record_human){
	delete $print_human{$i};
}

open HBV_OUT, ">$hbv_out" or die $!;
open HUMAN_OUT, ">$human_out" or die $!;

print HBV_OUT "reads_id\t\t\t\t\tref\tcigar\tpos\n";
print HUMAN_OUT "reads_id\t\t\t\t\tref\tcigar\tpos\n";

for my $key(keys %print_human){
	my @a = (split /\s+/, $print_human{$key})[0..3];
    my $p = join "\t", @a;
#	if($key =~ /_/){																	## modify at 11:49 2012-01-10
#		print HUMAN_OUT "repeat_reads_$p\n";											## modify at 11:49 2012-01-10
#	}																					## modify at 11:49 2012-01-10		
	print HUMAN_OUT "$p\n";
}	

for my $key(keys %print_hbv){
	print HBV_OUT "$print_hbv{$key}\n";
}

close HBV_OUT;
close HUMAN_OUT;
