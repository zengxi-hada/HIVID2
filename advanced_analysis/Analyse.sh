exec 1>output.debug
exec 2>output.someError
echo $1,$2
perl /public/home/chshen/auto_Image_Pro/get_all_txt.pl $1 $2;
/public/home/chshen/auto_Image_Pro/Run_pre.sh
/public/home/chshen/auto_Image_Pro/getHBV_CIRCOS.sh
/public/home/chshen/auto_Image_Pro/getNTCIRCOS.sh
/public/home/chshen/auto_Image_Pro/run_cpg_tfbs_tss.sh


perl /public/home/chshen/auto_Image_Pro/get_region.pl normal.txt tumor.txt

perl /public/home/chshen/auto_Image_Pro/check_bknum_bychr_partsample_noMXY_noN_telomere.pl normal.txt /public/home/chshen/auto_Image_Pro/hg19.info >normal.lib.chrNum
perl /public/home/chshen/auto_Image_Pro/check_bknum_bychr_partsample_noMXY_noN_telomere.pl tumor.txt /public/home/chshen/auto_Image_Pro/hg19.info > tumor.lib.chrNum
perl /public/home/chshen/auto_Image_Pro/result.pl normal.lib.chrNum>normal_result
perl /public/home/chshen/auto_Image_Pro/result.pl tumor.lib.chrNum>tumor_result
perl /public/home/chshen/auto_Image_Pro/getchr.pl
rm -f normal.lib.chrNum tumor.lib.chrNum normal_result tumor_result 


perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/Common_fragile_sites.xls.merge.wl normal.txt Normal_Common.distri
perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/NFRS.txt normal.txt Normal_NFRS.distri
perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/Rare_fragile  normal.txt Normal_Rare.distri
perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/Common_fragile_sites.xls.merge.wl tumor.txt Tumor_Common.distri
perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/NFRS.txt tumor.txt Tumor_NFRS.distri
perl /public/home/chshen/auto_Image_Pro/counk_bk_dis2.pl /public/home/chshen/auto_Image_Pro/Rare_fragile  tumor.txt Tumor_Rare.distri

perl /public/home/chshen/auto_Image_Pro/getRunFra.pl

rm -f Normal_NFRS.distri Normal_Common.distri Normal_Rare.distri Tumor_Common.distri Tumor_NFRS.distri Tumor_Rare.distri 

rm -f  autoRunFra.R autoRun.R config_NT_hmgene_g_gte10_2.conf2now hg19_CpGIsland.txt.len.pct2.500 NormalGene.circos NormalGene.circos.p normal.lib.chrNum normal_result normal.txt normal_virus.txt nownosub.conf RunCPG.R TumorGene.circos tumor.lib.chrNum tumor_result tumor.txt tumor_virus.txt
