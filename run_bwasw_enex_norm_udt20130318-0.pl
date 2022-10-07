#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
my $bin=dirname (abs_path ($0));

=head1 function
  
  This program is designed to autorun BWA-MEM and to grab breakpoint from the result of BWA-MEM.

=head1 usage
  
  perl $0
		-fq1	<str>		abstract path of human-un assembly fq file (paired-end after trimmomatic)
		-fq2	<str>		abstract path of hbv-un assembly fq file (paired-end after trimmomatic)
		-fq3	<str>		abstract path of un-un assembly fq file	(paired-end after trimmomatic)
		-fq4	<str>		abstract path of se-se assembly fq file (paired-end after trimmomatic)
		-fq5	<str>       abstract path of human-un assembly fq file (single-end after trimmomatic)
		-fq6	<str>       abstract path of hbv-un assembly fq file (single-end after trimmomatic)
		-fq7	<str>       abstract path of un-un assembly fq file (single-end after trimmomatic)
		-fq8	<str>       abstract path of se-se assembly fq file (single-end after trimmomatic)
		-fa1	<str>		abstract path of human ref []
		-fa2 	<str>		abstract path of HBV ref []
		-t	<int>		threshold value of alignment length in bwa mem [default 30]
		-o	<str>		abstract path of output directory
		-qsub			qsub job automatically or not, default not. it's must be used with "-vf"
		-vf	<num>		the memory need to be setted when qsub job automatically
		-filter	<str>		filter the un-unique seq
                -l	<str>		name of the sample

=head1 author 
  
  zengxi

=head1 version

  version1.0: 2011-11-10
  version1.2: 2012-04-06
  version2.0: 2021-02-18

=head1 last update
	2020-12

=cut	

my $fq1;
my $fq2;
my $fq3;
my $fq4;
my ($fq5, $fq6, $fq7, $fq8);
my ($qsub, $vf, $filter,$sample_name,$len);
my $fa1 = "/ifs1/ST_REHEAL/USER/PUBLIC_database/database/Homo_sapiens/HG19_noRandom/bwa_index/human.fa";
#my $fa2 = "/ifs1/ST_REHEAL/USER/zengxi/HBV/bwa/hbv_bwa_index/HBV.fa";
my $fa2 = "/ifs1/ST_REHEAL/USER/zengxi/save/HPV/mask_repeat/hpv_new.fa";
my $outdir;
my $threshold = 20;

GetOptions(
	'fq1=s' => \$fq1,
	'fq2=s' => \$fq2,
	'fq3=s' => \$fq3,
	'fq4=s' => \$fq4,
	'fq5=s' => \$fq5,
        'fq6=s' => \$fq6,
        'fq7=s' => \$fq7,
        'fq8=s' => \$fq8,
	'fa1=s' => \$fa1,
	'fa2=s' => \$fa2,
	't=i' => \$threshold,
	'qsub' => \$qsub,
	'vf=f' => \$vf,
	'o=s' => \$outdir,
	'filter' => \$filter,
	'l=s' => \$sample_name,
	's=i' => \$len,
);	


die `pod2text $0` unless ($fq1 && $fq2 && $fq3 && $fq4 && $fq5 && $fq6 && $fq7 && $fq8 && $outdir);
my $bwa = "$bin/bwa";
my $select_mutual_reads = "$bin/select_mutual_reads_enex.pl";
my $search_bk_human = "$bin/search_breakpoint_forHuman.pl";
my $search_bk_hbv = "$bin/search_breakpoint_forHBV_20130318.pl";
my $merge_bk_forHuman = "$bin/merge_breakpoint.pl";
my $merge_bk_forHBV = "$bin/merge_breakpoint_for_HBV_norm_20130318.pl";
my $merge_bk_forHuman_sm = "$bin/merge_breakpoint_sm.pl";
my $merge_bk_forHBV_sm = "$bin/merge_breakpoint_for_HBV_sm_norm_20130318.pl";
my $rm_dup = "$bin/rm_duplication.pl";
my $replacement = "$bin/check_replacement.pl";
my $hbv_len = "$bin/check_HBV_len.pl";
my $samtools ="samtools";
my $rm_assemblydup = "$bin/rm_assemblydup.pl";
my $msort = "$bin/msort";
my $alu_distance = "$bin/alu_distribution.pl";
my $annotation = "$bin/annotation1_addintergenic.pl";
my $gff = "$bin/Human_hg19_Refgene.gff";
my $alu_region = "$bin/nestedRepeats.txt.alu.gz";
my $circos = "$bin/circos";
my $process_pairwise = "$bin/process_pairwise_24s.pl";
my $circos_input_forRH = "$bin/circos_input_forRH.pl";
my $RHsup_human = "$bin/config_RHsup_human.conf";
my $gather_window = "$bin/gather_windows2_for_Human.pl";
my $circos_input_for_cluster = "$bin/circos_input_for_cluster.pl";
my $hg19_info = "$bin/hg19.info";
my $rm_pcr_dup = "$bin/restrict_rm_dup_only_consider_human-align_HZAU.pl";
my $check_chimera = "$bin/check_chimera.pl";
my $rm_pcr_dup_again = "$bin/try-fast.v3.py3.py";  #writen by zhouyi

mkdir $outdir if (not -e $outdir);
my $basename;
$basename = basename($fq1);
#my $basename2 = basename($fq2);
##if ($base =~ /_L\d+_(.*)_\d\.fq/){$basename = $1;}
$basename =~ s/reads_assembly_mergefa_quality_//;
$basename =~ s/\.clean_R1\.fq\.trimmo\.paired//;
$basename =~ s/\.clean_1\.fq\.trimmo\.paired//;
$basename =~ s/\.paired//;
#$basename = $base;

mkdir "$outdir/fq" if (not -e "$outdir/fq/");

#my $fq1_base = basename($fq1);

my $dir_fq1 = dirname($fq1);
my $base_fq1 = basename($fq1);
my $dir_fq2 = dirname($fq2);
my $base_fq2 = basename($fq2);
my $dir_fq3 = dirname($fq3);
my $base_fq3 = basename($fq3);
my $dir_fq4 = dirname($fq4);
my $base_fq4 = basename($fq4);

$base_fq1 =~ s/reads_assembly_mergefa_quality_//;
$base_fq2 =~ s/reads_assembly_mergefa_quality_//;
$base_fq3 =~ s/reads_assembly_mergefa_quality_//;
$base_fq4 =~ s/reads_assembly_mergefa_quality_//;

my ($left_fq1_1, $left_fq1_2, $left_fq2_1, $left_fq2_2, $left_fq3_1, $left_fq3_2, $left_fq4_1, $left_fq4_2);
$left_fq1_1 = "$dir_fq1/left_$base_fq1";
$left_fq1_2 = $left_fq1_1; $left_fq1_2 =~ s/1\.fq/2\.fq/; $left_fq1_2 =~ s/1\.fastq/2\.fastq/;

$left_fq2_1 = "$dir_fq2/left_$base_fq2";
$left_fq2_2 = $left_fq2_1; $left_fq2_2 =~ s/1\.fq/2\.fq/; $left_fq2_2 =~ s/1\.fastq/2\.fastq/;

$left_fq3_1 = "$dir_fq3/left_$base_fq3";
$left_fq3_2 = $left_fq3_1; $left_fq3_2 =~ s/1\.fq/2\.fq/; $left_fq3_2 =~ s/1\.fastq/2\.fastq/;

$left_fq4_1 = "$dir_fq4/left_$base_fq4";
$left_fq4_2 = $left_fq4_1; $left_fq4_2 =~ s/1\.fq/2\.fq/; $left_fq4_2 =~ s/1\.fastq/2\.fastq/;

open OUT, ">$outdir/bwa_mem_and_call_integration_sites.sh" or die $!;

print OUT "#!/bin/bash
#PBS -N HIVID_stp4
#PBS -l nodes=1:ppn=5
#PBS –l walltime=100:00:00
#PBS –l mem=12G
#PBS -q batch
#PBS -V
cd \$PBS_O_WORKDIR\n\n";

print OUT "date >&2\n";

print OUT "echo \"##  gunzip the fq file and merge them\" >&2\n";


mkdir "$outdir/human" unless -e "$outdir/human";
mkdir "$outdir/virus" unless -e "$outdir/virus";

print OUT "perl -p -i.bak -w -e 's/^\\\@trim_pe/\\\@left_trim_pe/g' $left_fq1_1 $left_fq1_2 $left_fq2_1 $left_fq2_2 $left_fq3_1 $left_fq3_2 $left_fq4_1 $left_fq4_2\n";     ### mark the reads not assembled successfully
print OUT "cat $fq1 $left_fq1_1 $left_fq1_2 $fq2 $left_fq2_1 $left_fq2_2 $fq3 $left_fq3_1 $left_fq3_2 $fq4 $left_fq4_1 $left_fq4_2 $fq5 $fq6 $fq7 > $outdir/fq/$basename.fq\n";  ###  modify at 22:02 2019-06-05

print OUT "echo \"##  gzip the merged fq file and remove duplication again\" >&2\n";  ### modify at 15:21 2012-02-19
print OUT "gzip -f $outdir/fq/$basename.fq\n";
print OUT "perl $rm_dup -a1 $outdir/fq/$basename.fq.gz -o $outdir/fq\n";      ### modify at 15:21 2012-02-19

print OUT "echo \"##  align the treated fq file with BWA-MEM\" >&2\n";
#print OUT "$bwa index $outdir/fq/index/$basename.fq\n";    ###  modify at 11:09 2011-11-11
print OUT "$bwa mem $fa1 $outdir/fq/rmdup_$basename.fq.gz > $outdir/human/human_$sample_name.sam\n";        ### modify at 15:21 2012-02-19
print OUT "$bwa mem $fa2 $outdir/fq/rmdup_$basename.fq.gz > $outdir/virus/virus_$sample_name.sam\n";            
print OUT "$samtools view -h -q 9 -F 4 -F 256 $outdir/human/human_$sample_name.sam > $outdir/human/human_$sample_name.uniq_map.sam\n";   ### filter the mutiple mapped reads to reserve the unique mapping reads

mkdir "$outdir/human/breakpoint" unless -e "$outdir/human/breakpoint";
mkdir "$outdir/virus/breakpoint" unless -e "$outdir/virus/breakpoint";

print OUT "echo \"##  search precise breakpoint\" >&2\n";
my ($preposition_dir ,$library_name) = $outdir=~/^(.*)\/step4\/(.*)$/;                                               ### modify at 17:03 2012-02-19
my @stat = <$preposition_dir/step3/$library_name/*.stat>; 															  ### modify at 17:03 2012-02-19
my @uniq_rate_file_arr = <$preposition_dir/step3/$library_name/SOAP/*.uniq_rate>;
my $uniq_rate_file = $uniq_rate_file_arr[0]; 

print OUT "perl $select_mutual_reads $outdir/virus/virus_$sample_name.sam $outdir/human/human_$sample_name.uniq_map.sam $outdir/virus/breakpoint/$sample_name\_mutual_sam.virus.endo $outdir/human/breakpoint/$sample_name\_mutual_sam.human.endo $outdir/virus/breakpoint/$sample_name\_mutual_sam.virus $outdir/human/breakpoint/$sample_name\_mutual_sam.human\n";

my $fl = "";
if($filter){$fl = "y"}else{$fl = "n";}

print OUT "perl $search_bk_hbv $outdir/virus/breakpoint/$sample_name\_mutual_sam.virus $outdir/virus/breakpoint/initialized_bk_$sample_name.virus $threshold $fl\n";
print OUT "perl $search_bk_human $outdir/human/breakpoint/$sample_name\_mutual_sam.human $outdir/human/breakpoint/initialized_bk_$sample_name.human $threshold $fl\n";

# check reads which may be rearrangement (嵌合/chimera) in human genome, rather than virus integration
print OUT "perl $check_chimera $outdir/human/human_$sample_name.sam $outdir/virus/virus_$sample_name.sam $outdir/human/breakpoint/initialized_bk_$sample_name.human $outdir/virus/breakpoint/initialized_bk_$sample_name.virus $outdir/human/breakpoint/use.initialized_bk_$sample_name.human $outdir/virus/breakpoint/use.initialized_bk_$sample_name.virus $outdir/human/breakpoint/low_confident.initialized_bk_$sample_name.human $outdir/virus/breakpoint/low_confident.initialized_bk_$sample_name.virus\n\n";

## process the high confident reads
print OUT "echo \"## process the high confident reads\" >&2\n";
print OUT "perl $merge_bk_forHBV $outdir/virus/breakpoint/use.initialized_bk_$sample_name.virus $stat[0] $outdir/virus/breakpoint/$sample_name\_virus_bk.final.stp1 $outdir/virus/breakpoint/$library_name.eff.stp y n $uniq_rate_file\n";    ### modify at 2019-07-04
print OUT "perl $merge_bk_forHuman $outdir/human/breakpoint/use.initialized_bk_$sample_name.human $stat[0] $outdir/human/breakpoint/$sample_name\_human_bk.final.stp1 $outdir/human/breakpoint/$library_name.eff.stp y n $uniq_rate_file\n";   ### modify at 2019-07-04
print OUT "$msort -k 1 -k n2 $outdir/virus/breakpoint/$sample_name\_virus_bk.final.stp1 > $outdir/virus/breakpoint/$sample_name\_virus_bk.final.stp1.sort; $msort -k 1 -k n2 $outdir/human/breakpoint/$sample_name\_human_bk.final.stp1 > $outdir/human/breakpoint/$sample_name\_human_bk.final.stp1.sort\n";    ### modify at 16:48 2012-04-06
print OUT "perl $merge_bk_forHBV_sm $outdir/virus/breakpoint/$sample_name\_virus_bk.final.stp1.sort $stat[0] $outdir/virus/breakpoint/$sample_name\_virus_bk.final.stp2 $outdir/virus/breakpoint/$library_name.eff.stp y 20\n";   ### modify at 14:00 2013-03-18
print OUT "perl $merge_bk_forHuman_sm $outdir/human/breakpoint/$sample_name\_human_bk.final.stp1.sort $stat[0] $outdir/human/breakpoint/$sample_name\_human_bk.final.stp2 $outdir/human/breakpoint/$library_name.eff.stp y 20\n";    ### modify at 14:16 2012-04-06
#print OUT "perl $rm_assemblydup $outdir/fq/rmdup_$sample_name.fq.gz $outdir/virus/breakpoint/$sample_name\_virus_bk.final.stp2 $outdir/virus/breakpoint/$sample_name\_virus_bk.final.stp2.uniq\n"; ### modify at 14:01 2013-03-18
print OUT "perl $rm_pcr_dup $outdir/human/breakpoint/$sample_name\_human_bk.final.stp2 $outdir/virus/breakpoint/$sample_name\_virus_bk.final.stp2 $outdir/human/human_$sample_name.uniq_map.sam $outdir/virus/virus_$sample_name.sam $outdir/human/breakpoint/$sample_name\_human_bk.final.stp2.uniq $outdir/virus/breakpoint/high_confident_$sample_name\_virus_bk.final.stp2.uniq\n";
## remove PCR duplicate again
#print OUT "python $bin/try-fast.py -fq1 /public/home/yzhou/test-uniq/infect-fastq/\$i.1.fq.gz -fq2 /public/home/yzhou/test-uniq/infect-fastq/\$i.2.fq.gz -i $outdir/human/breakpoint/$sample_name\_human_bk.final.stp2.uniq -o $outdir/human/breakpoint/$sample_name\_human_bk.final.stp2.uniq2 -id $sample_name\n";
(my $sample_name2 = $sample_name)=~s/3\./4\./;
($sample_name2 = $sample_name2)=~s/1\./2\./;
($sample_name2 = $sample_name2)=~s/1_/2_/;
##print $sample_name." ".$sample_name2."\n";
#print $sample_name2."\n";
#print OUT "python $bin/Uniq2.py -fq1 $outdir/../../step2/$sample_name/$sample_name.paired1.gz -fq2 $outdir/../../step2/$sample_name/$sample_name2.paired1.gz -i $outdir/human/breakpoint/$sample_name\_human_bk.final.stp2.uniq -o $outdir/human/breakpoint/$sample_name\_human_bk.final.stp2.uniq2 -id $sample_name -ref $bin/ref.list\n";

print OUT "gzip -d $outdir/../../step3/$sample_name/station_pair_end/$sample_name\_se_se.gz\n";
print OUT "perl -lane 'if (\$F[0]=~/1\$/) {print \$_.\"\\t+\\t$len\";} else {print \$_.\"\\t-\\t$len\";}' $outdir/../../step3/$sample_name/station_pair_end/$sample_name\_se_se >$outdir/../../step3/$sample_name/station_pair_end/$sample_name\_se_se_tran\n";
print OUT "mv -f $outdir/../../step3/$sample_name/station_pair_end/$sample_name\_se_se_tran $outdir/../../step3/$sample_name/station_pair_end/$sample_name\_se_se\n";
print OUT "gzip -f $outdir/../../step3/$sample_name/station_pair_end/$sample_name\_se_se\n";

## process the low confident reads which may be rearrangement rather than virus integration (the result of this part should not be used)
print OUT "echo \"## process the low confident reads which may be rearrangement rather than virus integration\" >&2\n";
print OUT "perl $merge_bk_forHBV $outdir/virus/breakpoint/low_confident.initialized_bk_$sample_name.virus $stat[0] $outdir/virus/breakpoint/low_confident.$sample_name\_virus_bk.final.stp1 $outdir/virus/breakpoint/low_confident.$library_name.eff.stp y n $uniq_rate_file\n";    ### modify at 2019-07-04
print OUT "perl $merge_bk_forHuman $outdir/human/breakpoint/low_confident.initialized_bk_$sample_name.human $stat[0] $outdir/human/breakpoint/low_confident.$sample_name\_human_bk.final.stp1 $outdir/human/breakpoint/low_confident.$library_name.eff.stp y n $uniq_rate_file\n";   ### modify at 2019-07-04
print OUT "$msort -k 1 -k n2 $outdir/virus/breakpoint/low_confident.$sample_name\_virus_bk.final.stp1 > $outdir/virus/breakpoint/low_confident.$sample_name\_virus_bk.final.stp1.sort; $msort -k 1 -k n2 $outdir/human/breakpoint/low_confident.$sample_name\_human_bk.final.stp1 > $outdir/human/breakpoint/low_confident.$sample_name\_human_bk.final.stp1.sort\n";    ### modify at 16:48 2012-04-06
print OUT "perl $merge_bk_forHBV_sm $outdir/virus/breakpoint/low_confident.$sample_name\_virus_bk.final.stp1.sort $stat[0] $outdir/virus/breakpoint/low_confident.$sample_name\_virus_bk.final.stp2 $outdir/virus/breakpoint/low_confident.$library_name.eff.stp y 20\n";   ### modify at 14:00 2013-03-18
print OUT "perl $merge_bk_forHuman_sm $outdir/human/breakpoint/low_confident.$sample_name\_human_bk.final.stp1.sort $stat[0] $outdir/human/breakpoint/low_confident.$sample_name\_human_bk.final.stp2 $outdir/human/breakpoint/low_confident.$library_name.eff.stp y 20\n";    ### modify at 14:16 2012-04-06
print OUT "perl $rm_pcr_dup $outdir/human/breakpoint/low_confident.$sample_name\_human_bk.final.stp2 $outdir/virus/breakpoint/low_confident.$sample_name\_virus_bk.final.stp2 $outdir/human/human_$sample_name.uniq_map.sam $outdir/virus/virus_$sample_name.sam $outdir/human/breakpoint/low_confident.$sample_name\_human_bk.final.stp2.uniq $outdir/virus/breakpoint/low_confident.$sample_name\_virus_bk.final.stp2.uniq\n\n";

#print OUT "$samtools view -b -h -S $outdir/virus/virus_$sample_name.sam -o $outdir/virus/virus_$sample_name.bam; $samtools view -b -h -S $outdir/human/human_$sample_name.uniq_map.sam -o $outdir/human/human_$sample_name.uniq_map.bam\n";
#print OUT "rm -f $outdir/virus/virus_$sample_name.sam $outdir/human/human_$sample_name.sam $outdir/human/human_$sample_name.uniq_map.sam $outdir/fq/$sample_name\_human_un.fq $outdir/fq/$sample_name\_hbv_un.fq $outdir/fq/$sample_name\_un_un.fq $outdir/fq/$sample_name\_se_se.fq $outdir/fq/$sample_name.fq.gz $outdir/fq/sort_$sample_name.fq\n";
#print OUT "gzip -f $dir_fq1/* $dir_fq2/* $dir_fq3/* $dir_fq4/*\n\n";

print OUT "date >&2\n";
#print OUT "echo \"All work has done!\" >&2\n";

if($qsub){
	$vf=$vf."g";
	system "qsub -cwd -l vf=$vf -o $outdir -e $outdir $outdir/bwa_mem_$sample_name.sh";
}
close OUT;
