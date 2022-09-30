perl /public/home/chshen/auto_Image_Pro/check_bknum_bychr_partsample_noMXY_noN_telomere.pl normal.txt hg19.info >normal.lib.chrNum
perl /public/home/chshen/auto_Image_Pro/check_bknum_bychr_partsample_noMXY_noN_telomere.pl tumor.txt hg19.info > tumor.lib.chrNum
perl /public/home/chshen/auto_Image_Pro/result.pl normal.lib.chrNum>normal_result
perl /public/home/chshen/auto_Image_Pro/result.pl tumor.lib.chrNum>tumor_result
perl /public/home/chshen/auto_Image_Pro/getchr.pl
