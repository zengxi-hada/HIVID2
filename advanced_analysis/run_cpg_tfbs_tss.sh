perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/hg19_CpGIsland.txt tumor.txt tumor.bk.cpg.distr
perl /public/home/chshen/auto_Image_Pro/counk_dis_win.pl -i tumor.bk.cpg.distr -w 50 -o tumor_normal.bk.0.cpg.distr.num.win50
head -21 tumor_normal.bk.0.cpg.distr.num.win50 > tumor_normal.bk.0.cpg.distr.num.win50.1000
perl -lane 'my $a=$F[-1]/10; print "$F[0]\t$F[1]\t$a"' tumor_normal.bk.0.cpg.distr.num.win50.1000 > tumor_normal.bk.0.cpg.distr.num.win50.1000.ts
head -11 tumor_normal.bk.0.cpg.distr.num.win50.1000.ts > tumor_normal.bk.0.cpg.distr.num.win50.500


perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/hg19_CpGIsland.txt normal.txt normal.bk.cpg.distr
perl /public/home/chshen/auto_Image_Pro/counk_dis_win.pl -i  normal.bk.cpg.distr -w 50 -o tumor_normal.bk.1.cpg.distr.num.win50
head -21 tumor_normal.bk.1.cpg.distr.num.win50 > tumor_normal.bk.1.cpg.distr.num.win50.1000
perl -lane 'my $a=$F[-1]/10; print "$F[0]\t$F[1]\t$a"' tumor_normal.bk.1.cpg.distr.num.win50.1000 > tumor_normal.bk.1.cpg.distr.num.win50.1000.ts
head -11 tumor_normal.bk.1.cpg.distr.num.win50.1000.ts > tumor_normal.bk.1.cpg.distr.num.win50.500
head -11 /public/home/chshen/auto_Image_Pro/hg19_CpGIsland.txt.len.pct2.ts > hg19_CpGIsland.txt.len.pct2.500


perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/tfbsConsSites.txt tumor.txt  tumor.bk.tfbs.distr
perl /public/home/chshen/auto_Image_Pro/counk_dis_win.pl -i tumor.bk.tfbs.distr -w 50 -o tumor_normal.bk.0.tfbs.distr.num.win50
head -21 tumor_normal.bk.0.tfbs.distr.num.win50 > tumor_normal.bk.0.tfbs.distr.num.win50.1000
perl -lane 'my $a=$F[-1]/10; print "$F[0]\t$F[1]\t$a"' tumor_normal.bk.0.tfbs.distr.num.win50.1000 > tumor_normal.bk.0.tfbs.distr.num.win50.1000.ts
head -11 tumor_normal.bk.0.tfbs.distr.num.win50.1000.ts>tumor_normal.bk.0.tfbs.distr.num.win50.500

perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/tfbsConsSites.txt normal.txt normal.bk.tfbs.distr
perl /public/home/chshen/auto_Image_Pro/counk_dis_win.pl -i normal.bk.tfbs.distr -w 50 -o tumor_normal.bk.1.tfbs.distr.num.win50
head -21 tumor_normal.bk.1.tfbs.distr.num.win50> tumor_normal.bk.1.tfbs.distr.num.win50.1000
perl -lane 'my $a=$F[-1]/10; print "$F[0]\t$F[1]\t$a"' tumor_normal.bk.1.tfbs.distr.num.win50.1000>tumor_normal.bk.1.tfbs.distr.num.win50.1000.ts
head -11 tumor_normal.bk.1.tfbs.distr.num.win50.1000.ts>tumor_normal.bk.1.tfbs.distr.num.win50.500




perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/switchDbTss.txt tumor.txt tumor.bk.tss.distr
perl /public/home/chshen/auto_Image_Pro/counk_dis_win.pl -i tumor.bk.tss.distr  -w 50 -o tumor_normal.bk.0.tss.distr.num.win50
head -21 tumor_normal.bk.0.tss.distr.num.win50>tumor_normal.bk.0.tss.distr.num.win50.1000
perl -lane 'my $a=$F[-1]/10; print "$F[0]\t$F[1]\t$a"' tumor_normal.bk.0.tss.distr.num.win50.1000>tumor_normal.bk.0.tss.distr.num.win50.1000.ts
head -11 tumor_normal.bk.0.tss.distr.num.win50.1000.ts >tumor_normal.bk.0.tss.distr.num.win50.500

perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/switchDbTss.txt normal.txt normal.bk.tss.distr
perl /public/home/chshen/auto_Image_Pro/counk_dis_win.pl -i normal.bk.tss.distr -w 50 -o tumor_normal.bk.1.tss.distr.num.win50
head -21 tumor_normal.bk.1.tss.distr.num.win50 > tumor_normal.bk.1.tss.distr.num.win50.1000
perl -lane 'my $a=$F[-1]/10; print "$F[0]\t$F[1]\t$a"' tumor_normal.bk.1.tss.distr.num.win50.1000>tumor_normal.bk.1.tss.distr.num.win50.1000.ts
head -11 tumor_normal.bk.1.tss.distr.num.win50.1000.ts > tumor_normal.bk.1.tss.distr.num.win50.500

output=`Rscript /public/home/chshen/auto_Image_Pro/merge3-2_modify.R`
perl /public/home/chshen/auto_Image_Pro/getPic.pl
rm -f ./*distr*
