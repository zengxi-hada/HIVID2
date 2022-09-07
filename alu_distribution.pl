#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;


my $alu_region = shift;
my $bk = shift;

die "perl $0 <alu_region> <bk>\n" unless ($alu_region && $bk);
if(-B $alu_region){
	open ALU, "<:gzip", $alu_region or die $!;
}else{
	open ALU, $alu_region or die $!;
}

if(-B $bk){
    open BK, "<:gzip", $bk or die $!;
}else{
    open BK, $bk or die $!;
}

my (%start_pos, %alu);

while(<ALU>){
	next if (!/Alu/i);
	chomp;
	my @a = split;
	my $chr = $a[1]; 
	my $start = $a[2];
	my $end = $a[3];	
	push @{$start_pos{$chr}}, $start;
	$alu{$chr}{$start} = $end;
}
close ALU;

#my $count = @start_pos;

while(<BK>){
	chomp;
	next if /NNSS/;
	my @a = split;
	my $chr = $a[2];
	my $pos = $a[3];
	$pos =~ s/,//g;
	my $rreads = pop @a;
	my $lreads = pop @a;
##	pop @a;
##	my $rreads = "";
##	my $lreads = "";
	my $count = @{$start_pos{$chr}};
	my $print_module1 = join ("\t", @a);
	my @search_array = sort {$a<=>$b} @{$start_pos{$chr}};
	my ($leftf, $rightf) = &BinarySearch(\@search_array, $pos);
	my $leftf_mate = $alu{$chr}{$start_pos{$chr}->[$leftf]};
	if($leftf == -1){                            ### the value is smaller than the first item of the array (-1, 0)
		my $distance = $start_pos{$chr}->[0] - $pos;
		print "$print_module1\t$distance\t$chr\_$start_pos{$chr}->[$leftf]\_$start_pos{$chr}->[$rightf]\t$lreads\t$rreads\n";
#		print "$distance\n";
		next;
	}elsif($leftf == $count-1){                  ### the value is larger than the last item of the array (count-1, $count)
		my $distance;
		if($leftf_mate >= $pos){
			my $distance = 0;
		}else{
			$distance = $pos - $leftf_mate;
		}
		print "$print_module1\t$distance\t$chr\_$start_pos{$chr}->[$leftf]\_$start_pos{$chr}->[$rightf]\t$lreads\t$rreads\n";
#		print "$distance\n";
		next;
	}
	if($leftf_mate>=$pos || $leftf==-100){          
		print "$print_module1\t0\t$chr\_$start_pos{$chr}->[$leftf]\_$start_pos{$chr}->[$rightf]\t$lreads\t$rreads\n";
#		print "0\n";
	}elsif($leftf_mate < $pos){
		my $distance = &min(abs($pos-$leftf_mate), abs($start_pos{$chr}->[$rightf]-$pos));
		print "$print_module1\t$distance\t$chr\_$start_pos{$chr}->[$leftf]\_$start_pos{$chr}->[$rightf]\t$lreads\t$rreads\n";
#		print "$distance\n";
	}
}

################################ subrontine ################################

sub BinarySearch  ########### retrun left and right offset ##########
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
                        return -100, $middle;					#### hit and -100 is a marker with no other sense
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


sub min
{
    my @a = @_;
    my $min = $a[0];
    for my $i(@a){
        if($i < $min){
            $min = $i;
        }
    }
    return $min;
}	
