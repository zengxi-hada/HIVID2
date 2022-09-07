#!usr/bin/perl -w
use strict;
use File::Basename;
use Cwd 'abs_path';
my $usage="perl $0 <list> <output file>";
die $usage if @ARGV!=2;
open OUT,">$ARGV[1]" or die $!;
print OUT "Sample name\tFC\tLane\tLibrary\tUseful_length\tInsertsize\tTotalReads\tTotalBases\tQ20\tGC%\tAdapter rate\tLow quality rate\tDuplication rate\tEffect reads\tHBV alignment\tHuman alignment\n";
open (IN,$ARGV[0]) or die "can't open $ARGV[0]:\n$usage\n";
while (<IN>){
	chomp;
	my @a=split /\s+/;
	open REPORT,"$a[-1]" or die "can't open $a[-1]\n";
	my ($TotalReads,$TotalBases,$Q20,$GC);
	while (my $report=<REPORT>){
		chomp $report;
		next if $report=~/#/;
		my @b=split /\s+/,$report;
		$TotalReads=$b[-1] if $b[-2]=~/TotalReads/;
		$TotalBases=$b[-1] if $b[-2]=~/TotalBases/;
		$Q20=$b[-1] if $b[-2]=~/Q20/;
		$GC=$b[-1] if $b[-2]=~/GC/;
	}
	close REPORT;
	my ($adapter_num,$lowqual_num);
	open RMA,"$a[-2]" or die "can't open $a[-2]\n";
	while (my $rma=<RMA>){
		chomp $rma;
		my @c=split /:/,$rma;
		$adapter_num=$c[1] if $c[0]=~/Adapter reads number/;
		$adapter_num=~s/\s//g  if $c[0]=~/Adapter reads number/;
		$lowqual_num=$c[1] if $c[0]=~/Low-quality reads number/;
		$lowqual_num=~s/\s//g if $c[0]=~/Low-quality reads number/;
		
	}
	close RMA;
	my $duplication;
	open (RMD,$a[-3]) or die "can't open $a[-3]\n";
	while(my $rmd=<RMD>){
		chomp $rmd;
		if ($rmd=~/duplication reads number/){
			$duplication=(split /:/,$rmd)[1];
			$duplication=~s/\s+//g;
		}
	}
	close RMD;
	open (HBV,$a[-4]) or die "can't open $a[-4]\n";
	my ($HBV_align,$Human_align);
	while (my $hbv=<HBV>){
		chomp $hbv;
		$HBV_align=(split /\t/,$hbv)[1] if $hbv=~/HBV percent/;
		$Human_align=(split /\t/,$hbv)[1] if $hbv=~/Human percent/;
	}
	close HBV;
	my $allread=$TotalReads*2;
	$TotalReads=sprintf("%.2f",$TotalReads/1000000);
	if($TotalReads>=1000){
		$TotalReads=sprintf("%.2f",$TotalReads/1000);
		$TotalReads="$TotalReads"."G";
	}else{
		 $TotalReads="$TotalReads"."M";	
	}
	$TotalBases=sprintf("%.2f",$TotalBases/1000000);
	if($TotalBases>=1000){
                $TotalBases=sprintf("%.2f",$TotalBases/1000);
                $TotalBases="$TotalBases"."G";
        }else{
                 $TotalBases="$TotalBases"."M";
        }
	print "$allread\n";
	print OUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$TotalReads;$TotalReads\t$TotalBases\t$Q20\t$GC\t";
	print OUT "$adapter_num(";
	printf OUT "%.2f",$adapter_num/$allread*100;
	print OUT "%)\t";
	print OUT "$lowqual_num(";
	printf OUT "%.2f",$lowqual_num/$allread*100;
	print OUT "%)\t";
	print OUT "$duplication(";
	printf OUT "%.2f",$duplication/$allread*100;
	print OUT "%)\t";
	print OUT $allread-$adapter_num-$lowqual_num-$duplication,"(";
	printf OUT "%.2f",($allread-$adapter_num-$lowqual_num-$duplication)/$allread*100;
	print OUT "%)\t";
	print OUT  "$HBV_align\t$Human_align\n";
}
close OUT;
close IN;
