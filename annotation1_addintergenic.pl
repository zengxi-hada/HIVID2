#!/usr/bin/perl

use strict;
use warnings;

my $gff = shift;
my $bk = shift;

die "perl $0 <.gff> <.bk>\n" unless ($gff && $bk);

open GFF, $gff or die $!;

my %record_gff;
my %record_ingene;
my ($id, $name);

while(<GFF>){
	chomp;
	my @a = split;
	my $chr = $a[0];
	my $spos = $a[3];
	my $epos = $a[4];
	my $style = $a[2];
	if(/mRNA/){
		$id = $a[-2];
		$id =~ s/;//;
		$id =~ s/ID=//;
		$name = $a[-1];	
		$name =~ s/;//;
		$name =~ s/name=//;
		push @{$record_gff{$chr}{"$spos\t$epos"}}, "$name:$id";
	}else{
		push @{$record_ingene{$chr}{"$spos\t$epos"}}, "$name:$id:$style";
	}
}
close GFF;

open BK, $bk or die $!;
while(<BK>){
    chomp;
    my @a = split;
	my $sample = $a[0];
    my $chr = $a[2];
    my $pos = $a[3];
	my $flag1 = 0; my $flag2 = 0;
	print "$sample\t$chr\t$pos\t";
	for my $i(keys %record_gff){
		next if $chr ne $i;
		for my $k(keys %{$record_gff{$i}}){
			my @value = @{$record_gff{$i}{$k}};
			my @value2 = map {"$_:Promoter"} @value;
			my $gene = join (";", @value2);
			my $spos = (split /\t/, $k)[0];
			my $epos = (split /\t/, $k)[1];
			if($pos <= $epos+1000 && $pos > $epos){
				my @value2 = map {"$_:Downstream"} @value;
		        my $gene = join (";", @value2);
				print "$gene;";
				$flag1 = 1;
			}elsif($pos >= $spos-5000 && $pos < $spos){
				my @value2 = map {"$_:Promoter"} @value;
				my $gene = join (";", @value2);
				print "$gene;";
				$flag1 = 1;
			}
		}
	}
	for my $i(keys %record_ingene){
		next if $chr ne $i;
		for my $k(keys %{$record_ingene{$i}}){
			my @value = @{$record_ingene{$i}{$k}};
			my $gene = join (";", @value);
			my $spos = (split /\t/, $k)[0];
            my $epos = (split /\t/, $k)[1];
			if($pos < $epos && $pos >= $spos){
                print "$gene;";
				$flag2 = 1;
            }
		}
	}
	if($flag1==0 && $flag2==0){
		my @tmp_keys = keys %{$record_gff{$chr}};
		my @search_array = map {[split /\t/, $_]->[0]} @tmp_keys;
		my %tmp_hash = map {[split /\t/, $_]->[0], [split /\t/, $_]->[1]} @tmp_keys;
		@search_array = sort {$a<=>$b} @search_array;
		my $count = @search_array;
		my ($left_offset, $right_offset) = &BinarySearch(\@search_array, $pos);
		my $left_offset_value = ($left_offset>=0) ? $search_array[$left_offset] : "nogene";
		my $right_offset_value = ($right_offset<$count) ? $search_array[$right_offset] : "nogene";
		my $left_offset_rightmate = ($left_offset>=0) ? $tmp_hash{$search_array[$left_offset]} : "nogene";
		my $right_offset_rightmate = ($right_offset<$count) ? $tmp_hash{$search_array[$right_offset]} : "nogene";
#		print "count of \@search_array: $count\n";
#		print "\n\$left_offset\t$left_offset\t\$right_offset\t$right_offset\n";
#		print "\$left_offset_value\t$left_offset_value\t\$right_offset_value\t$right_offset_value\n";
#		print "\$left_offset_rightmate\t$left_offset_rightmate\t\$right_offset_rightmate\t$right_offset_rightmate\n";
		my @sort_chr_array = sort {[split /\t/, $a]->[0] <=> [split /\t/, $b]->[0]} keys %{$record_gff{$chr}};
#		print "secondary key of \%record_gff (left): $left_offset_value\t$left_offset_rightmate\n";
#		print "secondary key of \%record_gff (right): $right_offset_value\t$right_offset_rightmate\n";
##		my $up_gene = ($left_offset>=0) ? join ("/", @{$record_gff{$chr}{"$left_offset_value\t$left_offset_rightmate"}}) : $sort_chr_array[0];
		my $up_gene = ($left_offset>=0) ? join ("/", @{$record_gff{$chr}{"$left_offset_value\t$left_offset_rightmate"}}) : "none";
#		my $down_gene = ($right_offset<$count) ? join ("/", @{$record_gff{$chr}{"$right_offset_value\t$right_offset_rightmate"}}) : $sort_chr_array[$#sort_chr_array];
		my $down_gene = ($right_offset<$count) ? join ("/", @{$record_gff{$chr}{"$right_offset_value\t$right_offset_rightmate"}}) : "none";
		my $up_distance = ($left_offset>=0) ? ($pos - $left_offset_rightmate) : "none";
		my $down_distance = ($right_offset<$count) ? ($right_offset_value - $pos) : "none";
		print "intergenic;  $up_gene:$up_distance;  $down_gene:$down_distance";
	}
	print "\n";
}
close BK;


##################################### sub routine ################################

sub BinarySearch
{
        my ($array1, $value) = @_;
        my @array = @$array1;
        my $middle = 0;
        my $left = 0;
        my $right = @array - 1;
        while($left <= $right)
        {
                $middle = ($left + $right) / 2;
                $middle = int($middle);
                if($array[$middle] == $value)
                {
                        return $middle, $middle;
                }
                elsif($array[$middle] < $value)
                {
                        $left = $middle + 1;
                }
                else
                {
                        $right = $middle - 1;
                }
        }
        return $right, $left;
}
