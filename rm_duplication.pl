#!usr/bin/perl -w
use strict;
use Getopt::Long;
use PerlIO::gzip;
use File::Basename;
my $usage=<<USAGE;
	PE perl $0 -a1 < fq1 file> -a2 < fq2 file> -o <outdir> -save -S <sort RAM)
	SE perl $0 -a1 < fq file> -o <outdir> -save 
	-save: whether save the duplication read 
USAGE
my($fq1_file,$fq2_file,$out,$save,$help,$ram);
GetOptions (
	'a1=s'=>\$fq1_file,
	'a2=s'=>\$fq2_file,
	'o=s'=>\$out,
	'S=s'=>\$ram,
	'save!'=>\$save,
	'help|?'=>\$help,
);
die $usage if (!defined $fq1_file|!defined $out);
$ram=2500000000 unless $ram;
die $usage if ($help);
`mkdir $out` if !-d $out;
$ram=($ram=~/(\d+)k$/)?($1*10**3):
         ($ram=~/(\d+)m$/)?($1*10**6):
         ($ram=~/(\d+)g$/)?($1*10**9):
         ($ram=~/(\d+)$/)?($1):(die "Wrong max memory input!");
my $name=basename $fq1_file;
$name=~s/_\d.fq.gz\b//;
open STAT,">$out/rmdup_$name.stat" or die $!;
print STAT "Remove duplication.. START TIME: ",`date`;
#---------------------------------------------
my ($sum_num,$dup_num,$seqL)=(0,0);
if ($fq2_file){
	&rm_PE_duplication($fq1_file,$fq2_file,\$sum_num,\$dup_num,\$seqL);
}else{
	&rm_SE_duplication($fq1_file,\$sum_num,\$dup_num,\$seqL);
}

#----------------------------------------------
print STAT "Original reads number:\t",$sum_num,"\n";
print STAT "Original bases number:\t",$sum_num*$seqL,"\n";
print STAT "duplication reads number:\t",$dup_num,"\n";
print STAT "duplication reads rate:\t",($dup_num/$sum_num)*100,"\n";
print STAT "Modified bases number:\t",($sum_num-$dup_num)*$seqL,"\n";
print STAT "Remove duplication.. END TIME: ",`date`;
close STAT;
#------------remove SE duplication------------
sub rm_SE_duplication{
	my($file,$num1,$num2,$len)=@_;
	my ($rm_dup,$sort_file,$dup);
	my $file_na= basename $file;
	$rm_dup="$out/rmdup_$file_na";
	$dup="$out/dup_$file_na";
	$sort_file="$out/sort_$file_na";	
	$sort_file=~s/\.gz$// if $sort_file=~/\.gz$/;
	open MERGE,">$sort_file" or die $!;
	if ($file=~/.gz\b/){
		open IN,"gzip -dc $file|" or die $!;
	}else{
		open IN,$file or die $!;
	}
	$/="\n";
	while(<IN>){
		chomp;
		my ($read_id)=(split)[0]; 
		my $seq=<IN>;
		<IN>;
		my $qual=<IN>;
		chomp($seq,$qual);
		print MERGE "$read_id\t$seq\t$qual\n";
	}
	close IN;
	close MERGE;
	`sort -s -S $ram -k 2 $sort_file -o $sort_file`;
	$rm_dup="$rm_dup.gz" unless $rm_dup=~/.gz\b/;
	$dup ="$dup.gz"  unless $dup=~/.gz\b/;
	open OUT,"|gzip > $rm_dup" or die $!;
	open DUP,"|gzip > $dup" or die $! if $save;
	my ($ID,$seq,$qual)=(0,0,0);
	open INT,"$sort_file" or die $!;
	$/="\n";
	while(<INT>){
		chomp;
		$$num1++;
		my @b=split;
		$$len=length $b[1];
		if ($seq eq $b[1]){
			print DUP "$b[0]\n$b[1]\n+\n$b[2]\n" if ($save);
			$$num2++;
		}else{
			print OUT "$b[0]\n$b[1]\n+\n$b[2]\n";
			$ID=$b[0];
                        $seq=$b[1];
                        $qual=$b[2];
		}
	}
#	print OUT "$ID\n$seq\n+\n$qual\n";
	close INT;
	close DUP if $save;
	close OUT;
##	`rm -f $sort_file`;
}
#--------- remove PE duplication -------------------------
sub rm_PE_duplication{
	my ($file1,$file2,$num1,$num2,$len)=@_;
	my ($rm_dup1,$dup1,$rm_dup2,$dup2,$sort_file);
	$file1=~/.*\//;
	my $name1 = $';
	$name1 = "$name.gz" unless $name1=~/\.gz\b/;
	$rm_dup1="$out/rmdup_$name1";
	$dup1="$out/dup_$name1";
	$file2=~/.*\//;
	my $name2=$';
	$name2="$name2.gz" unless $name2=~/\.gz\b/;
	$rm_dup2="$out/rmdup_$name2";
	$dup2="$out/dup_$name2";
	$name2=~s/_\d\w+\b//;
	$sort_file="$out/sort_$name2";
	$sort_file=~s/_\d.\w+.gz\b//;
	if($file1=~/.gz\b/){
		open FQ1,"gzip -dc $file1|" or die $!;
	}else{
		open FQ1,"$file1" or die $!;
	}
	if($file2=~/.gz\b/){
                open FQ2,"gzip -dc $file2|" or die $!;
        }else{
                open FQ2,"$file2" or die $!;
        }
	open MERGE," > $sort_file" or die $!;
	$/="@";
	<FQ1>;
	<FQ2>;
	while (<FQ1>){
		chomp;
		chomp(my $fq2=<FQ2>);
		my @a=split /\n/;
#		print $a[0],"\n";
		my @b=split /\n/,$fq2;
		print MERGE "$a[0]\t$b[0]\t$a[1]\t$b[1]\t$a[3]\t$b[3]\n";
	}
	close FQ1;
	close FQ2;
	`sort -s -S $ram -k 3,4 $sort_file -o $sort_file`;
	my ($ID,$SEQ,$QUAL,$id,$seq,$qual)=(0,0,0,0,0,0);
	open OUT1,"|gzip > $rm_dup1" or die $!;
	open OUT2,"|gzip > $rm_dup2" or die $!;
	open DUP1,"|gzip > $dup1" or die $! if $save;
	open DUP2,"|gzip > $dup2" or die $! if $save;
	open INT,"$sort_file" or die $!;
	$/="\n";
	while (<INT>){
		chomp;
		$$num1++;
		my @c=split;
		$$len=length $c[3];
		if ($SEQ eq $c[2] && $seq eq $c[3]){
			print DUP1 "@"."$c[0]\n$c[2]\n+\n$c[4]\n" if ($save);
			print DUP2 "@"."$c[1]\n$c[3]\n+\n$c[5]\n" if ($save);
			$$num2++;
		}else{
			print OUT1 "@"."$c[0]\n$c[2]\n+\n$c[4]\n";
                        print OUT2 "@"."$c[1]\n$c[3]\n+\n$c[5]\n";
			$ID=$c[0];
                        $id=$c[1];
                        $SEQ=$c[2];
                        $seq=$c[3];
                        $QUAL=$c[4];
                        $qual=$c[5];
		}
	}
	close OUT1;
	close OUT2;
	close DUP1 if $save;
	close DUP2 if $save;
	$$num1=$$num1*2;
	$$num2=$$num2*2;
	`rm -rf $sort_file`;
}
