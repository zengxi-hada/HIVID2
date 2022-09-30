cd Normal_pre
perl /public/home/chshen/auto_Image_Pro/sort.pl
perl /public/home/chshen/auto_Image_Pro/annotation1_addintergenic_liren_10k_new.pl /public/home/xzeng/bin/BGI_bin/annotation/Human_hg19_Refgene.gff break_all.txt >normal.anno
perl /public/home/chshen/auto_Image_Pro/sort2.pl normal.anno
mv all_bk_sort.anno normal.txt
cd ../

cd Tumor_pre
perl /public/home/chshen/auto_Image_Pro/sort.pl
perl /public/home/chshen/auto_Image_Pro/annotation1_addintergenic_liren_10k_new.pl /public/home/xzeng/bin/BGI_bin/annotation/Human_hg19_Refgene.gff break_all.txt >tumor.anno
perl /public/home/chshen/auto_Image_Pro/sort2.pl tumor.anno
mv all_bk_sort.anno tumor.txt
cd ../

cd Normal_pre_virus
perl /public/home/chshen/auto_Image_Pro/sort_virus.pl
mv break_all.txt normal_virus.txt
cd ../

cd Tumor_pre_virus
perl /public/home/chshen/auto_Image_Pro/sort_virus.pl
mv break_all.txt tumor_virus.txt
cd ../

cp ./Normal_pre/normal.txt ./Tumor_pre/tumor.txt ./Normal_pre_virus/normal_virus.txt ./Tumor_pre_virus/tumor_virus.txt ./
