#!usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
my $bin=dirname (abs_path ($0));

############################################################################################################
##	This program is to generate the folders and shell scripts for running step2	and step3				####
##  Author: Yi Shang 																					####
##			Zeng Xi																						####
##	Last update: 2020-12																				####
############################################################################################################

my $usage="perl $0 -o <out dir> -step <2 | 3> -list <FC Insertsize list> -c <config file> -f <the number of fqdata you want to find [10-37]> -filter <whether filter nounique align to human> -qsub -vf <RAM> -help\n";
my ($out,$step,$list,$con,$find,$qsub,$vf,$help,$filter);
GetOptions (
	'o=s'=>\$out,				# output directory
	'step=i'=>\$step,
	'list=s'=>\$list,
	'c=s'=>\$con,
	'f=s'=>\$find,
	'qsub!'=>\$qsub,
	'vf=s'=>\$vf,
	'help|?'=>\$help,
	'filter!'=>\$filter,
);

die $usage if $help;
die $usage if !$out|!$step|!$list|!$con;
if($vf){$vf=$vf."g" unless $vf=~/g/;}
`mkdir $out` if !-d $out;
die "can't mkdir $out\n" if !-d $out;
$find="20-37" if !$find;

#my $findfq="/ifs1/ST_REHEAL/USER/yishang/HBV/bin/find_fastq.pl";
my $findfq="$bin/find_fastq.pl";

#my $rmadapter="/ifs1/ST_REHEAL/USER/yishang/soft/rm_adapter_dup/rm_adaptor.pl";
#my $rmdup="/ifs1/ST_REHEAL/USER/yishang/soft/rm_adapter_dup/rm_duplication.pl";

#my $rmadapter_dup="/ifs1/ST_REHEAL/USER/yishang/soft/rm_adapter_dup/run_rmadapter_rmdup_sh.pl";
my $rmadapter_dup="$bin/run_rmadapter_rmdup_sh.pl";
#my $soap_and_station="/ifs1/ST_REHEAL/USER/yishang/HBV/bin/HBV_insertion.pl";
my $soap_and_station="$bin/HBV_insertion.pl";

my ($step1dir,$step2dir,$step3dir)=("$out/step1","$out/step2","$out/step3");
#print "$sample_list\n";
my $sample_list="$step1dir/sample.list";
#print "$sample_list\n";
if ($step==1){
	`mkdir $step1dir` if !-d $step1dir;
	$sample_list="$step1dir/sample.list";
	`perl $findfq -f $find -i $list -o $sample_list`;
	&readmestep1("$step1dir/step1_readme");
}elsif ($step==2){
	die "Please check step 1\n" if !-e $sample_list;
	`mkdir $step2dir` if !-d $step2dir; 
#	open OUTSH, ">$step2dir/trimmomatic.sh" or die $!;
	open FLT_LIST, ">$step2dir/rmadapter_lowquality.list" or die $!;
#	print OUTSH "#!/bin/bash\n#PBS -N trimmomatic.sh\n#PBS -l nodes=1:ppn=5\n#PBS –l walltime=100:00:00\n#PBS –l mem=10G\n#PBS -q batch\n#PBS -V\ncd \$PBS_O_WORKDIR\n";
#	print OUTSH "date\n";
	print FLT_LIST "Sample  FC      Lane    Library Useful_length   Insertsize      Pathway\n";
	open LIST, "$sample_list" or die "can't open $sample_list\n";
	<LIST>;
	while(my $LIST=<LIST>){
		chomp $LIST;
		my @a=split /\s+/,$LIST;
		for my $file($a[-1],$a[-2]){
			die "Please check step 1\n" if !-e $file;
		}
		my $sample_name = $a[0];
        mkdir "$step2dir/$sample_name";
		my $base_name1 = basename($a[-2]); my $base_name2 = basename($a[-1]);
        $base_name1 =~ s/\.gz$//; $base_name2 =~ s/\.gz$//;
		open OUTSH,">$step2dir/$sample_name/trimmomatic.sh" or die $!;
		print OUTSH "#!/bin/bash\n#PBS -N trimmomatic.sh\n#PBS -l nodes=1:ppn=5\n#PBS –l walltime=100:00:00\n#PBS –l mem=10G\n#PBS -q batch\n#PBS -V\ncd \$PBS_O_WORKDIR\n";
	    print OUTSH "date\n";
		print OUTSH "java -jar $bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $a[-2] $a[-1] $step2dir/$sample_name/$base_name1.trimmo.paired.gz $step2dir/$sample_name/$base_name1.trimmo.unpaired.gz $step2dir/$sample_name/$base_name2.trimmo.paired.gz $step2dir/$sample_name/$base_name2.trimmo.unpaired.gz ILLUMINACLIP:$bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60\n";
		print OUTSH "java -jar $bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $a[-2] $a[-1] $step2dir/$sample_name/$base_name1.trimmo.paired1.gz $step2dir/$sample_name/$base_name1.trimmo.unpaired1.gz $step2dir/$sample_name/$base_name2.trimmo.paired1.gz $step2dir/$sample_name/$base_name2.trimmo.unpaired1.gz ILLUMINACLIP:$bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0\n";		print OUTSH "date\n";
##		chdir $step2dir;
##		`qsub $step2dir/trimmomatic.sh`;
##		my $trimmomatic_log = `ls $step2dir/trimmomatic.sh.e*`;
##       chomp $trimmomatic_log;
		print FLT_LIST "$sample_name\t$a[4]\t$a[5]\ttrimmomatic.sh.e*\t$step2dir/$sample_name/$base_name1.trimmo.unpaired.gz\t$step2dir/$sample_name/$base_name2.trimmo.unpaired.gz\t$step2dir/$sample_name/$base_name1.trimmo.paired.gz\t$step2dir/$sample_name/$base_name2.trimmo.paired.gz\n";
	}
#		print OUTSH "date\n";
	close LIST;
	close OUTSH;
	close FLT_LIST;
	&readmestep2("$step2dir/step2_readme");
}elsif($step==3){
	my $rmadapter_lowquality_list="$step2dir/rmadapter_lowquality.list";
	die "Please check step 1 and step 2\n" if !-e $rmadapter_lowquality_list;
	open LIST,"$rmadapter_lowquality_list" or die "can't open $rmadapter_lowquality_list\n";
	<LIST>;
	while(my $LIST=<LIST>){
                chomp $LIST;
                my @a=split /\s+/,$LIST;
				die "Please check $rmadapter_lowquality_list \"Insertsize\" \n" unless  $a[-6]=~/\d+/;
                for my $file($a[-1],$a[-2],$a[-3]){
                        die "Please check step 1 and step 2\n" if !-e $file;
                }
        }
	close LIST;
	`mkdir $step3dir` if !-d $step3dir;
	($filter)?`perl $soap_and_station -o $step3dir -list $rmadapter_lowquality_list -c $con -soap -station -filter`:`perl $soap_and_station -o $step3dir -list $rmadapter_lowquality_list -c $con -soap -station` if !$qsub;
#	print "perl $soap_and_station -o $step3dir -list $rmadapter_lowquality_list -c $con -soap -station -filter\n";
	($filter)?`perl $soap_and_station -o $step3dir -list $rmadapter_lowquality_list -c $con -soap -station -filter -qsub -vf $vf`:`perl $soap_and_station -o $step3dir -list $rmadapter_lowquality_list -c $con -soap -station -qsub -vf $vf` if $qsub;
	&readmestep3("$step3dir/step3_readme");
}else{
	print "Please check option \"-step\"\n";
}
#sub sample_name{
#	my ($dir,$readstep)=@_;
#	opendir (DIR,$dir)or die "can't open $dir \n";
#	my @file=readdir DIR;
#	closedir DIR;
#	for (@file){
#		next if $_=~/[\.]{1,2}/;
#		next if !-d "$dir/$_";
#		if ($readstep==2){
#			&readmestep2("$dir/$_/step2_readme");
#		}elsif ($readstep==3){
#			&readmestep3("$dir/$_/step3_readme");
#		}
#	}
#}
		
sub readmestep1{
	my ($step1read)=@_;
	open OUT,">$step1read" or die "can't open $step1read\n";
	my $purpose="Finding the fq file  of date\n";
	my $options_findfq="
Options:
    -f STR  the number of fqdata you want to find [10-37]
    -i STR  the file info list [ID Library FC]
    -o STR  output file
    -h Help Information"
;
	print OUT "Purpose\t $purpose method:\nfirst\twe should set up a list,the first row is the name of sample;the second row is the library of sample; the third row is the FC number of sample.\nsecond\t run \"$findfq\"\n perl $findfq\n $options_findfq\nQuestion:\nwhy can't the step find the fq file  of date?\n first, you should chick the list of samples' information,the first row is the name of sample;the second row is the library of sample; the third row is the FC number of sample.\nSecondly,-f Range of the option -f needs changing\n" ;
	close OUT;
}
sub readmestep2{
	my($step2read)=@_;
        open OUT,">$step2read" or die "can't open $step2read\n";
	print OUT "Puopose\trm adapter and duplication,low quality reads;\nMethod:\nrun outdir/step2/samples' name/sample_*.sh\nthe mermory of running script is 20g\n when the script have been runned,the file would be created that is list of file's name about removing adapter and duplication\nQuestion:\nHow to do if the script can't be runned?\n you should check whether exist the list file that was created after running the step1\n"; 
	close OUT;
}
sub readmestep3{
	my($step3read)=@_;
	open OUT,">$step3read" or die "can't open $step3read\n";
	print OUT "Purpose\t alignment,station and PER assembly for the fq file after removing adapter and duplication\nThere are three file (HBV_soap.sh,Human_soap.sh,station.sh)be created very sample after run step 3,HBV_soap.sh and Human_soap.sh is script that the data aligned to HBV and human reference sequence by soap,the mermory is about 6g when running the Human_soap.sh or HBV_soap.sh;the station.sh is statistic analysis and PER,the mermory is about 12g;when you run the station.sh ,you must wait out the HBV_soap.sh and Human_soap.sh were over,the HBV_soap.sh and Human_soap.sh are runned  at the same time\n";
	close OUT;
}


	
