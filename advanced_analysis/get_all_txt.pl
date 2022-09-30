#!usr/bin/perl
#
use strict;
use warnings;
use File::Basename;

my $dir_tumor=shift;
my $dir_normal=shift;
die "$0 <tumor.stp4> <normal.stp4>" unless ($dir_tumor && $dir_normal);
chomp(my $now_dir=`pwd`);
print "$dir_tumor,$dir_normal\n";

`mkdir Tumor_pre Normal_pre`;
`mkdir Tumor_pre_virus Normal_pre_virus Image`;
#`cp sort.pl sort2.pl Normal_pre`;
#`cp sort.pl sort2.pl Tumor_pre`;
`ls $dir_tumor/*/human/breakpoint/*.final |xargs -i cp {} $now_dir/Tumor_pre`;
`ls $dir_normal/*/human/breakpoint/*.final |xargs -i cp {} $now_dir/Normal_pre`;
`ls $dir_tumor/*/virus/breakpoint/[^chimera]*.uniq |xargs -i cp {} $now_dir/Tumor_pre_virus`;
`ls $dir_normal/*/virus/breakpoint/[^chimera]*.uniq |xargs -i cp {} $now_dir/Normal_pre_virus`;

#`cp sort_virus.pl Normal_pre_virus`;
#`cp sort_virus.pl Tumor_pre_virus`;
#`cd ./Normal_pre`;
#`perl sort.pl`;
#`cd ../Tumor_pre`;
#`perl sort.pl`;
#my $file=
#`./run.sh`;
#print "111\n";



