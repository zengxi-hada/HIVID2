#!/usr/bin/perl

use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
my $bin=dirname (abs_path ($0));

=head1 function

  This program is designed to run the flow of virus integration project automatically. In details, it's to detect the integration breakpoints in condition that virus DNA inserts into human genome. The output contains information about position of breakpoints, support reads for both human and virus. etc.

=head1 usage

  perl $0
			-o	<str>		absolute path of output directory
			-l	<str>		file of list for samples
			-stp	<str>		the steps of the flow:1 for finding fqs; 2 for removing pollution of fq data; 3 for running soap, 
						doing statistics and assembly; 4 for running bwasw
			-c	<str>		absolute path of config file while running soap
			-f	<int>		the number of fqdata you want to find [10-37]
			-qsub			qsub or not, default not. it must be used with the parameter 'vf'
			-vf2	<num>		the memory setted for step2 when qsub jobs automatically
			-vf3	<num>		the memory setted for step3 when qsub jobs automatically
			-vf4	<num>		the memory setted for step4 when qsub jobs automatically
			-fa1	<str>		the absolute path of human reference when performing bwa-mem [hg19]
			-fa2	<str>		the absolute path of  pathogene reference when performing bwa-mem [virus]
			-filter	<str>		choose uniq alignment reads or not in soap

=head1 author
  
  zengxi@mail.hzau.edu.cn	
  zengxi@genomics.org.cn
  yishang@genomics.org.cn
  chenshengpei@genomics.org.cn
  15071254117@sina.cn (shenchenhang)
  yzhou.zoe@qq.com (zhouyi)

=head1 version

  version1.0: 2011-11-10
  version1.1: 2012-02-15
  version1.2: 2013-03-18
  version2.0: 2019-06-20
  version2.1: 2020-11-22
=cut

my ($outdir, $list, $stp, $qsub, $vf2, $vf3, $vf4, $filter);
my $fqnum = 36;
my $config = "$bin/ConfigHBV_19";
#my $prepare_data = "/ifs1/ST_REHEAL/USER/zengxi/bin/HBV/yishang_bin/HBV_foruse.pl";
#my $prepare_data = "/ifs1/ST_REHEAL/USER/yishang/HBV/bin/HBV.pl";
#my $prepare_data = "/ifs1/ST_REHEAL/PMO/SZY11044_HBV/bin/yishang_bin/HBV.pl";
my $prepare_data = "$bin/HBV.pl";
my $search_point = "$bin/run_bwasw_enex_norm_udt20130318-0.pl";
my $fa1 = "/ifs1/ST_REHEAL/USER/zengxi/save/cervical_carcinoma/mask_index/human/hg19_hpv_mask_final.fa";
my $fa2 = "/ifs1/ST_REHEAL/USER/zengxi/save/cervical_carcinoma/mask_index/hpv/hpv_mask.fa";

GetOptions(
	'o=s' => \$outdir,
	'l=s' => \$list,
	'stp=s' => \$stp,
	'c=s' => \$config,
	'qsub' => \$qsub,
	'f=s'=> \$fqnum,
	'vf2=f' => \$vf2,
	'vf3=f' => \$vf3,
	'vf4=f' => \$vf4,
	'fa1=s' => \$fa1,
	'fa2=s' => \$fa2,
	'filter' => \$filter,
);

die `pod2text $0` unless ($outdir && $list && $stp);

if(!$filter){
#	print "zengxi\n";
	if($stp == 1){
		system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum";	
	}elsif($stp == 2){
		if($qsub){
			system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum -qsub -vf $vf2";
		}else{
			system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum";
			print "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum\n";
		}	
	}elsif($stp == 3){
		if($qsub){
			system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum -qsub -vf $vf3";
		}else{
			system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum";
		}
	}elsif($stp == 4){
		mkdir "$outdir/step4" unless -e "$outdir/step4";
		open IN,$list or die $!;
		my @a = <IN>;
		for my $i(@a){
			chomp $i;
			my @b = split /\s+/, $i;
			my $sample_name = $b[0];
			my $lib = $b[1];
			my $fc = $b[2];
			mkdir "$outdir/step4/$sample_name" unless -e "$outdir/step4/$sample_name";
			my @fq1 = <$outdir/step3/$sample_name/reads_assemble/human_un/reads_assembly_mergefa_quality*>;
			my @fq2 = <$outdir/step3/$sample_name/reads_assemble/hbv_un/reads_assembly_mergefa_quality*>;
			my @fq3 = <$outdir/step3/$sample_name/reads_assemble/un_un/reads_assembly_mergefa_quality*>;
			my @fq4 = <$outdir/step3/$sample_name/reads_assemble/se_se/reads_assembly_mergefa_quality*>;
			my @fq5 = <$outdir/step3/$sample_name/reads_assemble_single-end/human_un/*unpaired*>;
                        my @fq6 = <$outdir/step3/$sample_name/reads_assemble_single-end/hbv_un/*unpaired*>;
                        my @fq7 = <$outdir/step3/$sample_name/reads_assemble_single-end/un_un/*unpaired*>;
                        my @fq8 = <$outdir/step3/$sample_name/reads_assemble_single-end/se_se/*unpaired*>;
#			print $fq1[0], "\n";
#
			if($qsub){
				system "perl $search_point -fq1 $fq1[0] -fq2 $fq2[0] -fq3 $fq3[0] -fq4 $fq4[0] -fq5 $fq5[0] -fq6 $fq6[0] -fq7 $fq7[0] -fq8 $fq8[0] -fa1 $fa1 -fa2 $fa2 -o $outdir/step4/$sample_name -qsub -vf $vf4";
			}else{
#				print "$search_point,$fq1[0],$fq2[0],$fq3[0],$fq4[0],$fq5[0],$fq6[0],$fq7[0],$fq8[0],$fa1,$fa2\n";
				system "perl $search_point -fq1 $fq1[0] -fq2 $fq2[0] -fq3 $fq3[0] -fq4 $fq4[0] -fq5 $fq5[0] -fq6 $fq6[0] -fq7 $fq7[0] -fq8 $fq8[0] -fa1 $fa1 -fa2 $fa2 -o $outdir/step4/$sample_name -l $sample_name";
			}
		}
	}
}else{
	if($stp == 1){
        system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum -filter";
    }elsif($stp == 2){
        if($qsub){
            system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum -filter -qsub -vf $vf2";
        }else{
            system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum -filter";
        }
    }elsif($stp == 3){
        if($qsub){
            system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum -filter -qsub -vf $vf3";
        }else{
            system "perl $prepare_data -o $outdir -list $list -step $stp -c $config -f $fqnum -filter";
        }
	}elsif($stp == 4){
		my $len;
		open LIST,"$outdir/step1/sample.list" or die $!;
		while(<LIST>) {next if /Insertsize/; my @a=split;$len=$a[5];last;}
	#	print $len."\n";
        mkdir "$outdir/step4" unless -e "$outdir/step4";
        open IN,$list or die $!;
#        my @a = <IN>;
#        for my $i(@a){
		while(<IN>){
			chomp;
			my @a = split;
#            chomp $i;
            my $sample_name = $a[0];
		#	print "2222222222222$sample_name\n";
            my $lib = $a[1];
            my $fc = $a[2];
            mkdir "$outdir/step4/$sample_name" unless -e "$outdir/step4/$sample_name";
            my @fq1 = <$outdir/step3/$sample_name/reads_assemble_pair-end/human_un/reads_assembly_mergefa_quality*>;
            my @fq2 = <$outdir/step3/$sample_name/reads_assemble_pair-end/hbv_un/reads_assembly_mergefa_quality*>;
            my @fq3 = <$outdir/step3/$sample_name/reads_assemble_pair-end/un_un/reads_assembly_mergefa_quality*>;
            my @fq4 = <$outdir/step3/$sample_name/reads_assemble_pair-end/se_se/reads_assembly_mergefa_quality*>;
	    my @fq5 = <$outdir/step3/$sample_name/reads_unassemble_single-end/human_un/*unpaired*>;
	    my @fq6 = <$outdir/step3/$sample_name/reads_unassemble_single-end/hbv_un/*unpaired*>;
	    my @fq7 = <$outdir/step3/$sample_name/reads_unassemble_single-end/un_un/*unpaired*>;
	    my @fq8 = <$outdir/step3/$sample_name/reads_unassemble_single-end/se_se/*unpaired*>;
            if($qsub){
                system "perl $search_point -fq1 $fq1[0] -fq2 $fq2[0] -fq3 $fq3[0] -fq4 $fq4[0] -fq5 $fq5[0] -fq6 $fq6[0] -fq7 $fq7[0] -fq8 $fq8[0] -fa1 $fa1 -fa2 $fa2 -o $outdir/step4/$sample_name -qsub -vf $vf4 -filter";
            }else{
			#	print "sample".$sample_name."\n";
#				print "perl $search_point -fq1 $fq1[0] -fq2 $fq2[0] -fq3 $fq3[0] -fq4 $fq4[0] -fq5 $fq5[0] -fq6 $fq6[0] -fq7 $fq7[0] -fq8 $fq8[0] -fa1 $fa1 -fa2 $fa2 -o $outdir/step4/$sample_name -filter";
                system "perl $search_point -fq1 $fq1[0] -fq2 $fq2[0] -fq3 $fq3[0] -fq4 $fq4[0] -fq5 $fq5[0] -fq6 $fq6[0] -fq7 $fq7[0] -fq8 $fq8[0] -fa1 $fa1 -fa2 $fa2 -o $outdir/step4/$sample_name -filter -l $sample_name -s $len";
##				print "perl $search_point -fq1 $fq1[0] -fq2 $fq2[0] -fq3 $fq3[0] -fq4 $fq4[0] -fq5 $fq5[0] -fq6 $fq6[0] -fq7 $fq7[0] -fq8 $fq8[0] -fa1 $fa1 -fa2 $fa2 -o $outdir/step4/$sample_name -filter -l $sample_name -s $len\n";
#				print "perl $search_point -fq1 $fq1[0] -fq2 $fq2[0] -fq3 $fq3[0] -fq4 $fq4[0] -fa1 $fa1 -fa2 $fa2 -o $outdir/step4/$sample_name -filter\n";
#				print "perl $search_point -fq1 $fq1[0] -fq2 $fq2[0] -fq3 $fq3[0] -fq4 $fq4[0] -fq5 $fq5[0] -fq6 $fq6[0] -fq7 $fq7[0] -fq8 $fq8[0] -fa1 $fa1 -fa2 $fa2 -o $outdir/step4/$sample_name -filter\n";
#				print "$outdir/step3/$sample_name/reads_assemble_pair-end/human_un/\n";
            }
        }
    }
}
close IN;
