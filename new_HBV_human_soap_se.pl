#!usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use PerlIO::gzip;

#####################################################################################################################################
### this program is to extract the reads which may contain virus DNA based on soap results of unpaired reads after trimmomatic ######
###	Autor: Zeng Xi																											   ######
#####################################################################################################################################

my $usage="perl $0 -hs <human soap SE> -hu <human soap unmap> -bs <HBV soap SE> -bu < HBV soap unmap> -se <out SE_SE> -sb<out HBV_UN> -sh <out human_UN> -un <out UN_UN> -stat <out stat file> -f1 <read1> -o <read file dir> -filter <whether filter non-unique human alignment>\n";
my ($hu_se,$hu_un,$hb_se,$hb_un,$se_se,$hbv_un,$human_un,$un_un,$stat,$read1,$out,$filter);
GetOptions (
        'hs=s'=>\$hu_se,						# single end read maps to human genome
        'hu=s'=>\$hu_un,						# single end read cannot map to human genome
        'bs=s'=>\$hb_se,						# single end read maps to virus genome
        'bu=s'=>\$hb_un,						# single end read cannot map tp virus genome
        'se=s'=>\$se_se,						# the single end read can both map on human and virus genome
        'sb=s'=>\$hbv_un,						# the single end read can map onto virus genome, but cannot map onto human genome
        'sh=s'=>\$human_un,						# the single end read can map onto human genome, but cannot map onto virus genome
        'un=s'=>\$un_un,						# the single end read can not mapped onto human genome, but cannot map onto both genome
        'stat=s'=>\$stat,
        'f1=s'=>\$read1,						# the unpaired reads after trimmomatic treatment
        'o=s'=>\$out,
		'filter!'=>\$filter,					# whether filter non-unique human alignment (unique means one read map onto only onto one position)
);

die $usage if (!$hu_se|!$hu_un|!$hb_se|!$hb_un|!$se_se|!$hbv_un|!$human_un|!$un_un|!$stat);

my (%huse,%huun,%hbse,%hbun);
($filter)?&read_pese($hu_se,\%huse,"human"):&read_pese($hu_se,\%huse,"0");
&read_pese($hb_se,\%hbse,"0");
&read_un($hu_un,\%huun);
&read_un($hb_un,\%hbun);

=h
for my $i(keys %hbun){
	print "$i\n";
}
=cut

my (%idpp,%idss,%idbu,%idhu,%iduu);
my ($se_se_num,$se_un_num)=(0,0,0);

open OUT2, "|gzip > $se_se" or die $!;
open OUT3, "|gzip > $hbv_un" or die $!;
open OUT4, "|gzip > $human_un" or die $!;
open OUT5, "|gzip > $un_un" or die $!;

for my $id (keys %hbse){
	my $hbsechr=(split /\s+/,$hbse{$id})[0];
	if(exists $huse{$id}){
		print OUT2 "$id"."\t",$hbse{$id},"\n";
		my $husechr = (split /\s+/,$huse{$id})[0];
		$idss{$id}="$husechr"."_$hbsechr";
#		delete $huse{$id};
		$se_se_num++;
	}elsif(exists $huun{$id}){
		print OUT3 "$id"."\t",$hbse{$id},"\n";
		$idbu{$id}="unmap_$hbsechr";
#		delete $huun{$id};
		$se_un_num++;
	}
#	delete $hbse{$id};
}

my ($un_se_num,$un_un_num)=(0,0,0);
for my $id(keys %hbun){
#	print "$id\n";
	if(exists $huse{$id}){
		print OUT4 "$id"."\t",$huse{$id},"\n";
		my $husechr=(split /\s+/,$huse{$id})[0];
		$idhu{$id}="$husechr"."_unmap";
#		delete $huse{$id};
		$un_se_num++;
	}elsif(exists $huun{$id}){
		print OUT5 "$id"."\t",$huun{$id},"\n";
		$un_un_num++;
		$iduu{$id}="unmap_unmap";
#		delete $huun{$id};
	}
#	delete $hbun{$id};
}

close OUT2;
close OUT3;
close OUT4;
close OUT5;

my $totail = $se_se_num + $se_un_num + $un_se_num + $un_un_num;
#print "$se_se_num $se_un_num $un_se_num $un_un_num $totail\n";
open STAT,">>$stat" or die $!;
if($totail ==0) {$totail=1000000000;}
print STAT "\n\nFor unpaired reads after trimmomatic\n";
print STAT "SE(human)\tUnmap(human)\n";
my @SE=($se_se_num,$se_un_num);
percent("SE(virus)",$totail,@SE);
my @UN=($un_se_num,$un_un_num);
percent("Unmap(virus)",$totail,@UN);
print STAT "virus percent\t";
print STAT $se_se_num+$se_un_num,"(";
printf STAT "%0.2f",($se_se_num+$se_un_num)/$totail*100;
print STAT "%)\n";
print STAT "effective Ratio\t";
printf STAT "%0.2f",($se_se_num)/$totail*100;
print STAT "%\n";
print STAT "Human percent\t",$se_se_num+$un_se_num,"(";
printf STAT "%0.2f",($se_se_num+$un_se_num)/$totail*100;
print STAT "%)\n";
if($totail==1000000000) {$totail=0;}
print STAT "Unpaired-Total\t$totail\n";							# how many single reads (not pairs)
close STAT;

die if (!$read1|!$out);
open FQ1,"gzip -cd $read1|" or die "can't open $read1";
`mkdir $out` if !-e $out && !-d $out;
die "can't mkdir $out " if !-e $out;
my $dir_sese="$out/se_se";
my $dir_hbvun="$out/hbv_un";
my $dir_humanun="$out/human_un";
my $dir_unun="$out/un_un";
`mkdir $dir_sese $dir_hbvun $dir_humanun $dir_unun` if !-d $dir_sese|!-d $dir_hbvun|!-d $dir_humanun|!-d $dir_unun;
my $name_read1=basename ($read1);
$name_read1=~s/.gz\b//;

open SS1,">$dir_sese/$name_read1" or die "can't open $dir_sese/$name_read1";
open HBVU1,">$dir_hbvun/$name_read1" or die "can't open $dir_hbvun/$name_read1";
open HUMAN1,">$dir_humanun/$name_read1" or die "can't open $dir_humanun/$name_read1";
open UNUN1,">$dir_unun/$name_read1" or die "can't open $dir_unun/$name_read1";

#$/="@";
#<FQ1>;
while (<FQ1>){
		my $id_fq1_str = $_;
		$id_fq1_str =~ s/^\@//;
	#	my $id_fq1 = (split /\s+/, $id_fq1_str)[0];
        my $id_fq1 = (split /\s+/, $id_fq1_str)[0];
		(my $id_fq11)=$id_fq1=~m/(\/[12]$)/;
		$id_fq1=~s/\/[12]$//;
       # print "$id_fq11\n";
      #  $id_fq1=~s/[12]$//;
        my $seq_fq1 = <FQ1>;
        my $line3_fq1 = <FQ1>;
        my $qual_fq1 = <FQ1>;
#        my @d=split;
		my $fq1_str = $id_fq1_str.$seq_fq1.$line3_fq1.$qual_fq1;
        if(exists $idss{$id_fq1}){
             	print SS1 "@"."trim_se#$idss{$id_fq1}".$fq1_str;
        }elsif(exists $idbu{$id_fq1}){
                print HBVU1 "@"."trim_se#$idbu{$id_fq1}".$fq1_str;
        }elsif(exists $idhu{$id_fq1}){
                print HUMAN1 "@"."trim_se#$idhu{$id_fq1}".$fq1_str;
        }elsif(exists $iduu{$id_fq1}){
                print UNUN1 "@"."trim_se#$iduu{$id_fq1}".$fq1_str;
        }
}
close SS1;
close HBVU1;
close HUMAN1;
close UNUN1;


############################# subroutine ############################
sub read_pese{
	my ($soap_file,$hash,$ter)=@_;
	if($soap_file=~/.gz$/){
		open IN,"gzip -cd $soap_file|" or die "can't open $soap_file\n";
	}else{
		open IN,"$soap_file" or die "can't open $soap_file\n";
	}
	while(<IN>){
		chomp;
		my ($soap_id,$unq,$read_order,$chr,$align_pos)=(split /\s+/)[0,3,4,7,8];
		if($unq>1 and $ter eq "human"){
			next;
		}
		$soap_id=~s/\/[12]$//;
		$$hash{$soap_id} = "$chr\t$align_pos";
	}
	close IN;
}

sub read_un{
	my ($soap_file,$hash)=@_;
	if($soap_file=~/.gz$/){
                open IN,"gzip -cd $soap_file|" or die "can't open $soap_file\n";
        }else{
                open IN,"$soap_file" or die "can't open $soap_file\n";
        }
#	my ($id,$soap_un)=(0,0);
	while(<IN>){
		chomp;
		$_=~s/^>//;
		$_=~s/\/[12]$//;
		my $id = $_;
		<IN>;						# skip the sequence line
		$$hash{$id}=1;
#		print "$id\n";
	}
	close IN;
}

sub percent{
        my ($head,$total,@data)=@_;
		my $percent;
        print STAT "$head\t";
        for(@data){
				if($total == 0){$percent = 0;}
                else{
                     $percent=100*$_/$total;
                }

                print STAT "$_(";
                printf STAT  "%0.2f",$percent;
                print STAT "%)\t";
        }
        print STAT "\n";
}
