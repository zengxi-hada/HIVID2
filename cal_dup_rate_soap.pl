#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;

my $soap_file = shift;
my $sample_list = shift;				# sample.list in step2

die "perl $0 <soap.file> <read.len>\n" unless ($soap_file && $sample_list);

my $read_len;
open SL, $sample_list or die $!;
my %readlen_sample;
while(<SL>){
	next if(/Insertsize/);
	my @a = split;
	my $sample_id = $a[0];
	$read_len = (split /;/, $a[4])[0];
	$readlen_sample{$sample_id} = $read_len;
}
close SL;
#print "$read_len\n";


my ($sample_id) = $soap_file=~/step3\/(\w+)\/SOAP/;

my %h;
#open SF, $soap_file or die $!;
open SF,"gzip -cd $soap_file|" or die "can't open $soap_file\n";
#open SF, $soap_file or die "can't open $soap_file\n";
while(<SF>){
	chomp;
	my @a = split /\s+/, $_;			# the first line
	my $read_id1 = $a[0];
	my $chr1 = $a[7];
	my $pos1 = $a[8];
	my $len1 = $a[5];
	$read_id1 =~ s/\/1//;
	chomp (my $line2 = <SF>);			# the second line
	my @b = split /\s+/, $line2;        
    my $read_id2 = $b[0];
    my $chr2 = $b[7];
    my $pos2 = $b[8];
    my $len2 = $b[5];
	$read_id2 =~ s/\/2//;
#	print "$read_id\t$chr\t$pos\t$len\n";
#	$h{$read_id} = "$chr\t$pos\t$len" if($len == $read_len);
	next if ($read_id1 ne $read_id2);
	$read_len = $readlen_sample{$sample_id};
	if($read_len==$len1 && $read_len==$len2){
#	if(141==$len1 && 141==$len2){
		push @{$h{"$chr1:$pos1:$len1\t$chr2:$pos2:$len2"}}, $read_id1;		# read_id1 = read_id2
	}
}
close SF;
my ($total_count, $dup_count, $uniq_count) = (0, 0, 0);
for my $pos_str(keys %h){
	$total_count += @{$h{$pos_str}};
	$dup_count += $total_count - 1;
	$uniq_count += 1;
}
my $uniq_reads_rate=0;
if($total_count != 0)
{
	 $uniq_reads_rate = $uniq_count / $total_count;
}
else {
 $uniq_reads_rate = 1;
}

print "uniq_reads_rate\t$uniq_reads_rate\n";
	
