#!/user/bin/perl

use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
my $bin=dirname (abs_path ($0));

=head1 function

  This program is to generate a all-in-one shell script for HIVID2 pipeline.

=head1 usage

  perl all_in_one.pl
 		-o 		<str>		absolute path of output directory
 		-ts		<str>		total sample list
 		-fa1		<str>		the absolute path of human reference when performing bwa-mem [human]
 		-fa2		<str>		the absolute path of pathogene reference when performing bwa-mem [virus]
		-bin		<str>		the absolute path of HIVID2 program (optional, default is the path of all_in_one.pl)
		-c		<str>		the absolute path of config file for running soap

=head1 authors

  Zeng,Xi: zengxi@mail.hzau.edu.cn, Xi.Zeng@childrens.harvard.edu
  Yi,Shang: yishang@genomics.org.cn
  Chen,Shengpei: chenshengpei@genomics.org.cn
  Shen,chenhang: 15071254117@sina.cn
  Zhou,Yi: yzhou.zoe@qq.com
  Wang,Yuyouye:

=head1 version

  version1.0: 2011-11-10
  version1.1: 2012-02-15
  version1.2: 2013-03-18
  version2.0: 2019-06-20
  version2.1: 2020-11-22
  version2.2: 2022-11-02

=cut

my ($outdir, $bin_dir, $total_list, $config, $fa1, $fa2);
GetOptions(
    'o=s' => \$outdir,
	'bin=s' => \$bin_dir,
    'ts=s' => \$total_list,
    'c=s' => \$config,
    'fa1=s' => \$fa1,
    'fa2=s' => \$fa2,
);

die `pod2text $0` unless ($outdir && $total_list && $bin_dir && $config && $fa1 && $fa2);

open TL, $total_list or die $!;
while(<TL>){
	next if /#/;
	chomp;
	my @a = split;
	my $sample_id = $a[0];
	mkdir "$outdir";
	mkdir "$outdir/$sample_id";
	mkdir "$outdir/$sample_id/step1";
	open SL, ">$outdir/$sample_id/step1/sample.list" or die $!;
	open LT, ">$outdir/$sample_id/list" or die $!;
	print SL "Sample\tFC\tLane\tLibrary\tUseful_length\tInsertsize\tPathway\n";
	print SL "$_\n";
	print LT "$sample_id\t$sample_id\t$sample_id\n";
	close SL;
	close LT;

	open OS, ">$outdir/$sample_id/$sample_id\_all_in_one.sh" or die $!;
	print OS "echo \"\n## step 2\" >&2\n";
	print OS "perl $bin_dir/main.pl -o $outdir/$sample_id -l $outdir/$sample_id/list  -c $config -stp 2\n";
	print OS "sh $outdir/$sample_id/step2/$sample_id/trimmomatic.sh\n";

	print OS "\necho \"\n## step 3\" >&2\n";
	print OS "perl $bin_dir/main.pl -o $outdir/$sample_id -l $outdir/$sample_id/list  -c $config -stp 3\n";
	print OS "sh $outdir/$sample_id/step3/$sample_id/Human_virus_soap.sh\n";
	print OS "sh $outdir/$sample_id/step3/$sample_id/station.sh\n";

	print OS "\necho \"\n## step 4\" >&2\n";
	print OS "perl $bin_dir/main.pl -o $outdir/$sample_id -l $outdir/$sample_id/list  -c $config -stp 4 -filter -fa1 $fa1 -fa2 $fa2\n";
	print OS "sh $outdir/$sample_id/step4/$sample_id/bwa_mem_and_call_integration_sites.sh\n";
        print OS "sh $bin_dir/run_update_breakpoints.sh $outdir/$sample_id $sample_id $bin_dir $fa1 $fa2";
	close OS;
}
close TL;
