#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $usage=<<'USAGE';
Options:
    -f STR  the number of fqdata you want to find [10-37]
    -i STR  the file info list [ID Library FC]
    -o STR  output file
    -h Help Information
USAGE

my $fqdata="15-38";
my ($in,$out,$help);
GetOptions(
	'f=s' => \$fqdata,
	'i=s' => \$in,
	'o=s' => \$out,
	'h|?' => sub{die "$usage\n";},
);

die "$usage\n" unless (defined($in) && defined($out));
open INFO,"<$in" or die "Can't open input file $in!\n";
open OUT,">$out" or die "Can't open output file $out!\n";

my ($start,$end)=(split /-/,$fqdata)[0,-1];
my (@pe,@se);
while(<INFO>){
	chomp;
	my @line=split;
	my $sample=shift @line;
	my $library=join "*",reverse(@line);
	my $bool=0;
	for(my $i=$start;$i<=$end;$i++){
		my $dir="/share/fqdata$i";
		next unless (-d $dir);
		chomp(my @find=`find $dir -name "*$library*_1.fq.gz"`);
		if(@find){
			foreach my $fq1(@find){
				my ($fq,$fqdir)=fileparse($fq1);
				my ($fc,$lane,$lib)=(split "_",$fq)[2,3,4];
				(my $fq2=$fq1)=~s/_1.fq.gz/_2.fq.gz/;
				(my $rep=$fq1)=~s/_1.fq.gz/.report/;
				if(!-f $rep){
					chomp(my @rep=`find $fqdir -name "*report"`);
					if(@rep==1){$rep=$rep[0];}else{($rep)=grep {$_ =~ /$lib/} @rep;}
				}
				unless(-f $rep){
					print OUT "$sample\tWrong!\tNo report!\n";
					$bool=1;
					next;
				}
				my ($length,$insert,$sd)=&rep($rep);
				if(-f $fq2){
					push @pe, "$sample\t$fc\t$lane\t$lib\t$length\t$insert\t$fq1\t$fq2\t$rep\n";
				}else{
					push @se, "$sample\t$lib\t$fc\t$lane\t$fq1\n";
				}
			}
			$bool=1;
		}
	}
	print OUT "$sample\tWrong!\tNo such file!\n" unless $bool;
}
close INFO;

if(@pe){
	print OUT "Sample\tFC\tLane\tLibrary\tUseful_length\tInsertsize\tPathway\n";
	print OUT @pe;
}
if(@se){
	print OUT "Sample\tLibrary\tFC\tLane\tPathway\n";
	print OUT @se;
}
close OUT;

sub rep{
	open REP,"<",$_[0];
	my ($length,$insert,$sd)=("unknown","unknown","unknown");
	while(<REP>){
		chomp;
		if(/Length\b/){$length=(split)[-1];}
		if(/InsertSize\b/){$insert=(split)[-1];}
		if(/InsertSizeSD\b/){$sd=(split)[-1];}
	}
	close REP;
	return ($length,$insert,$sd);
}
