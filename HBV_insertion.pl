#!usr/bin/perl -w

#############################################################################################################
###  This program is to generate the folders and shell scripts for running step3              			 ####
###  Author: Yi Shang                                                                                    ####
###          Zeng Xi                                                                                     ####
###  Last update: 2020-12                                                                                ####
#############################################################################################################

use strict;
use Getopt::Long;
use PerlIO::gzip;
use File::Path;
use File::Basename;
use Cwd qw/abs_path/;
my $bin=dirname (abs_path ($0));

my $usage=" perl $0 -o <out dir> -list <sample.list> -c <Config file> -filter <whether filter nounique align to human> -soap  -station  -qsub -vf -h";
my ($out,$list,$con,$soap,$station,$qsub,$vf,$help,$filter);
GetOptions (
	'o=s'=>\$out,							# step3 dir
	'list=s'=>\$list,						# rmadapter_lowquality.list
	'c=s'=>\$con,
	'soap!'=>\$soap,
	'station!'=>\$station,					# classify the reads according alignment results
	'qsub!'=>\$qsub,
	'vf=s'=>\$vf,
	'filter!'=>\$filter,
	'h|help!'=>\$help,
);
die $usage if $help;
die $usage if (!$out|!$list);
if ($vf){$vf="$vf"."g" unless $vf=~/g/}
$out=abs_path($out);
mkpath($out);
open (CON,$con)or die $!;
print $con;
$/=">";
my %config;
<CON>;
#print "$config";
#open FILE,">/public/home/xzeng/project/Dongfang/chshen/RealDataHIVID/HBV_in_PLC/hivid2/test.txt" or die $!;
while(my $line=<CON>){
	chomp $line;
	#chop $line;
	$line=~s/\n$//;
#	print FILE $line."  line\n";
        my ($type,$thing)=split /\=/,$line;
		print($type);
        if($type eq "virus_config"){
                my @array=split /\n/,$thing;
                foreach (@array){
                        my ($read_len,$co)=split /:/,$_;
			next if !$read_len;
                        $config{"virus_config"}{$read_len}=$co;
                }
                next;
        }
	if($type eq "Human_config"){
                my @array=split /\n/,$thing;
                foreach (@array){
                        my ($read_len,$co)=split /:/,$_;
                        next if !$read_len;
                        $config{"Human_config"}{$read_len}=$co;
                }
                next;
        }
        $config{$type}=$thing;
}
close CON;

my $station_pl="$bin/new_HBV_human_soap.pl";
my $station_pl_se = "$bin/new_HBV_human_soap_se.pl";
my $station_break="$bin/station_breakpiont.pl";
my $cal_dup_rate = "$bin/cal_dup_rate_soap.pl";
my $step1_folder = $out;
$step1_folder =~ s/step3/step1/;

open (LIST,$list) or die $!;
#print FILE $list."\n";
#
$/="\n";
#open LS,">$out/allsample.list" or die $!;
<LIST>;
while(<LIST>){
	chomp;
	my @a=split /\s+/;
	my $sam=abs_path("$out/$a[0]");
	mkpath($sam);
	my ($info,$virus_info,$Human_info);
	if ($soap){
		my $read_len=$a[1];
		my($low,$max);
		$config{"insert_sd"}=~/-(\d+)/;
		$low=$1;
		$config{"insert_sd"}=~/\+(\d+)/;
		$max=$&;
		$low=$a[2]-$low*1.5;
		$max=$a[2]+$max*1.5;
		my $sample="$sam/SOAP";									# sample is folder path of the SOAP files
		`mkdir $sample` if !-d $sample;
		my $virus_soap_pair="$sample/virus\_$a[0].pair.soap";
		my $virus_soap_single="$sample/virus\_$a[0].single.soap";
		my $virus_soap_unmap="$sample/virus\_$a[0].unmap.soap";
		my $Human_soap_pair="$sample/Human\_$a[0].pair.soap";
		my $Human_soap_single="$sample/Human\_$a[0].single.soap";
		my $Human_soap_unmap="$sample/Human\_$a[0].unmap.soap";
		
		my $virus_soap_se = "$sample/virus\_$a[0].se.soap";
		my $virus_soap_se_unmap = "$sample/virus\_$a[0].unmap.se.soap";
		my $Human_soap_se = "$sample/Human\_$a[0].se.soap";
		my $Human_soap_se_unmap = "$sample/Human\_$a[0].unmap.se.soap";
		open HUMAN,">$sam/Human_virus_soap.sh" or die $!;
		print HUMAN "#!/bin/bash\n#PBS -N Human_virus_soap.sh\n#PBS -l nodes=1:ppn=5\n#PBS –l walltime=100:00:00\n##PBS –l mem=10G\n#PBS -q batch\n#PBS -V\ncd \$PBS_O_WORKDIR\n";
		print HUMAN "date\n";
#		print FILE "Human:@a\n";
		my $merge_trimmo_unpaired = $a[-3]; $merge_trimmo_unpaired =~ s/_R2//; $merge_trimmo_unpaired =~ s/_2//; $merge_trimmo_unpaired =~ s/\.2//;$merge_trimmo_unpaired =~ s/2\./\./;$merge_trimmo_unpaired =~ s/\.gz$//;
#		print HUMAN "if [ ! -f $merge_trimmo_unpaired.gz ]; then\n";
#		print HUMAN "gunzip $a[-4] $a[-3]; ";
		$a[-4] =~ s/\.gz//; $a[-3] =~ s/\.gz//;
		print HUMAN "cat $a[-4].gz $a[-3].gz > $merge_trimmo_unpaired.gz;\n\n";
		print HUMAN $config{"soap"}," -a $a[-2] -b $a[-1] -D ",$config{"ref_human"}," -o $Human_soap_pair -2 $Human_soap_single -u $Human_soap_unmap -m $low -x $max -p 8 ",$config{"virus_config"}{$read_len},"\n\n";
#		print $config{"soap"}," -a $a[-2] -b $a[-1] -D ",$config{"ref_human"}," -o $Human_soap_pair -2 $Human_soap_single -u $Human_soap_unmap -m $low -x $max -p 8 ",$config{"virus_config"}{$read_len},"\n\n";
		print HUMAN $config{"soap"}," -a $merge_trimmo_unpaired.gz -D ",$config{"ref_human"}," -o $Human_soap_se -u $Human_soap_se_unmap -p 8 ", $config{"virus_config"}{$read_len},"\n";
#		close HUMAN;
#		open VIRUS,">$sam/virus_soap.sh" or die $!;
#		print VIRUS "if [!-f $merge_trimmo_unpaired.gz];then\n";
#        print VIRUS "gunzip $a[-4] $a[-3]; ";
		$a[-4] =~ s/\.gz//; $a[-3] =~ s/\.gz//;
#		print VIRUS "cat $a[-4] $a[-3] > $merge_trimmo_unpaired; gzip $merge_trimmo_unpaired; gzip $a[-4] $a[-3]\nfi\n\n";
#		print VIRUS "#!/bin/bash\n#PBS -N virus_soap.sh\n##PBS -l nodes=1:ppn=10\n##PBS –l walltime=100:00:00\n##PBS –l mem=10G\n##PBS -q batch\n##PBS -V\n#cd \$PBS_O_WORKDIR\n";
		print HUMAN $config{"soap"}," -a $a[-2] -b $a[-1] -D ",$config{"ref_virus"}," -o $virus_soap_pair -2 $virus_soap_single -u $virus_soap_unmap -m $low -x $max -p 8 ", $config{"Human_config"}{$read_len},"\n\n";
		print HUMAN $config{"soap"}," -a $merge_trimmo_unpaired.gz -D ",$config{"ref_virus"}," -o $virus_soap_se -u $virus_soap_se_unmap -p 8 ",$config{"Human_config"}{$read_len},"\n\n";
##		print HUMAN "gzip -f $virus_soap_pair $virus_soap_single $virus_soap_unmap $Human_soap_pair $Human_soap_single $Human_soap_unmap $Human_soap_se $Human_soap_se_unmap $virus_soap_se $virus_soap_se_unmap\n\n";
		print HUMAN "gzip -f $sample/*soap\n\n";
		print HUMAN "date\n";
		close HUMAN;

#		$virus_info=`qsub -cwd -l vf=$vf "$sam/virus_soap.sh" -o $sample -e $sample` if $qsub;
#		$Human_info=`qsub -cwd -l vf=$vf "$sam/Human_soap.sh" -o $sample -e $sample` if $qsub;

		if($station){
			my $station_pair_end="$sam/station_pair_end";
			my $station_single_end="$sam/station_single_end";
			`mkdir $station_pair_end` if !-d $station_pair_end;
			`mkdir $station_single_end` if !-d $station_single_end;
			my $pe_pe= "$station_pair_end/$a[0]\_pe_pe.gz";
			my $se_se="$station_pair_end/$a[0]\_se_se.gz";
			my $hbv_un="$station_pair_end/$a[0]\_virus_un.gz";
			my $human_un="$station_pair_end/$a[0]\_Human_un.gz";
			my $un_un="$station_pair_end/$a[0]\_un_un.gz";
			
			my $se_se_single_end="$station_single_end/$a[0]\_se_se.gz";
            my $hbv_un_single_end="$station_single_end/$a[0]\_virus_un.gz";
            my $human_un_single_end="$station_single_end/$a[0]\_Human_un.gz";
            my $un_un_single_end="$station_single_end/$a[0]\_un_un.gz";

			my $stat="$sam/$a[0].stat";
			my $out_assemble="$sam/reads_assemble_pair-end";					# paired end reads
			my $out_assemble_se = "$sam/reads_unassemble_single-end";				# single end reads
			my $breakpoint="$sam/$a[0]\_breakpoint.xls";
			open STAT,">$sam/station.sh" or die $!;
			print STAT "#!/bin/bash\n#PBS -N station.sh\n#PBS -l nodes=1:ppn=5\n#PBS –l walltime=100:00:00\n#PBS –l mem=10G\n#PBS -q batch\n#PBS -V\ncd \$PBS_O_WORKDIR\n\n";
			print STAT "date\n";
#			print STAT "gzip -d $virus_soap_pair $virus_soap_single $virus_soap_unmap $Human_soap_pair $Human_soap_single $Human_soap_unmap\n"; #test fail need
			($filter)?print STAT "perl $station_pl -reads_assembly $bin/overlap_pair_trim.new -margefa_qua $bin/margefa_qua.pl -hp $Human_soap_pair.gz -hs $Human_soap_single.gz -hu $Human_soap_unmap.gz -bp $virus_soap_pair.gz -bs $virus_soap_single.gz -bu $virus_soap_unmap.gz -pe $pe_pe -se $se_se -sb $hbv_un -sh $human_un -un $un_un  -stat $stat -f1 $a[-2] -f2 $a[-1] -o $out_assemble -filter\n\n" : print STAT "perl $station_pl -reads_assembly $bin/overlap_pair_trim.new -margefa_qua $bin/margefa_qua.pl -hp $Human_soap_pair.gz -hs $Human_soap_single.gz -hu $Human_soap_unmap.gz -bp $virus_soap_pair.gz -bs $virus_soap_single.gz -bu $virus_soap_unmap.gz -pe $pe_pe -se $se_se -sb $hbv_un -sh $human_un -un $un_un  -stat $stat -f1 $a[-2] -f2 $a[-1] -o $out_assemble \n\n";
			($filter)?print STAT "perl $station_pl_se -hs $Human_soap_se.gz -hu $Human_soap_se_unmap.gz -bs $virus_soap_se.gz -bu $virus_soap_se_unmap.gz -se $se_se_single_end -sb $hbv_un_single_end -sh $human_un_single_end -un $un_un_single_end -stat $stat -f1 $merge_trimmo_unpaired.gz -o $out_assemble_se -filter\n\n":print STAT "perl $station_pl_se -hs $Human_soap_se.gz -hu $Human_soap_se_unmap.gz -bs $virus_soap_se.gz -bu $virus_soap_se_unmap.gz -se $se_se_single_end -sb $hbv_un_single_end -sh $human_un_single_end -un $un_un_single_end -stat $stat -f1 $merge_trimmo_unpaired.gz -o $out_assemble_se \n\n";
			print STAT "perl $cal_dup_rate $Human_soap_pair.gz $step1_folder/sample.list > $Human_soap_pair.uniq_rate\n";
##			print STAT "perl $station_break -i $se_se -h $human_un -l $a[5] -o $breakpoint\n\n";
#			print STAT "gzip -f $virus_soap_pair $virus_soap_single $virus_soap_unmap $Human_soap_pair $Human_soap_single $Human_soap_unmap\n";
			print STAT "date\n";
			close STAT;
			if($qsub){
				$virus_info=~/\d+/;
				my $virus_id=$&;
				$Human_info=~/\d+/;
				my $Human_id=$&;
				`qsub -hold_jid $virus_id,$Human_id  -cwd -l vf=$vf "$sam/station.sh" -o $sample -e $sample`if $qsub;
			}
		}
	}
}
close LIST;
my $DATA="$bin/HBV_Dataproduction.pl";

open DA,">$out/Dataproduction.sh" or die $!;
print DA "perl $DATA $out/allsample.list Dataproduction.xls\n";
close DA;
