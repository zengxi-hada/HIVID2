perl /public/home/chshen/auto_Image_Pro/gather_windows2_for_greed_accumulation_nosub.pl -bk normal_virus.txt -l /public/home/chshen/auto_Image_Pro/hbv.len -w 100 -ol 1 -o tumor_normal.virus.bk.n.w100.nosub
perl /public/home/chshen/auto_Image_Pro/gather_windows2_for_greed_accumulation_nosub.pl -bk tumor_virus.txt -l /public/home/chshen/auto_Image_Pro/hbv.len -w 100 -ol 1 -o tumor_normal.virus.bk.t.w100.nosub
perl /public/home/chshen/auto_Image_Pro/circos_input_for_cluster.pl tumor_normal.virus.bk.n.w100.nosub >tumor_normal.virus.bk.n.w100.nosub.circos
perl /public/home/chshen/auto_Image_Pro/circos_input_for_cluster.pl tumor_normal.virus.bk.t.w100.nosub >tumor_normal.virus.bk.t.w100.nosub.circos
perl /public/home/chshen/auto_Image_Pro/gather_windows2_for_samplefrq_nosub.pl -bk normal_virus.txt -l /public/home/chshen/auto_Image_Pro/hbv.len -w 100 -ol 1 -o tumor_normal.virus.bk.n.sfrq.w100.nosub
perl /public/home/chshen/auto_Image_Pro/gather_windows2_for_samplefrq_nosub.pl -bk tumor_virus.txt -l /public/home/chshen/auto_Image_Pro/hbv.len -w 100 -ol 1 -o tumor_normal.virus.bk.t.sfrq.w100.nosub
perl /public/home/chshen/auto_Image_Pro/circos_input_for_cluster.pl tumor_normal.virus.bk.n.sfrq.w100.nosub > tumor_normal.virus.bk.n.sfrq.w100.nosub.circos
perl /public/home/chshen/auto_Image_Pro/circos_input_for_cluster.pl tumor_normal.virus.bk.t.sfrq.w100.nosub > tumor_normal.virus.bk.t.sfrq.w100.nosub.circos
cp /public/home/chshen/auto_Image_Pro/nownosub.conf ./
circos -conf nownosub.conf
rm -f tumor_normal.virus.bk.n.w100.nosub tumor_normal.virus.bk.t.w100.nosub tumor_normal.virus.bk.n.w100.nosub.circos tumor_normal.virus.bk.t.w100.nosub.circos tumor_normal.virus.bk.n.sfrq.w100.nosub tumor_normal.virus.bk.t.sfrq.w100.nosub tumor_normal.virus.bk.n.sfrq.w100.nosub.circos tumor_normal.virus.bk.t.sfrq.w100.nosub.circos nownosub.conf
