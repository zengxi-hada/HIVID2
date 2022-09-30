#$0 normal $1 tumor
perl -lane 'print "$F[3]"' normal.txt >genenormal.txt
perl -lane 'print "$F[3]"' tumor.txt >genetumor.txt
sed -i 's#,#\n#g' genenormal.txt
sed -i 's#,#\n#g' genetumor.txt
perl -lane 'print $_ if /\w+/' genenormal.txt >genenormal1.txt
perl -lane 'print $_ if /\w+/' genetumor.txt >genetumor1.txt
sort genenormal1.txt >genenormal2.txt
sort genetumor1.txt >genetumor2.txt
perl /public/home/chshen/auto_Image_Pro/DataPro.pl genenormal2.txt >geneNormal.txt
perl /public/home/chshen/auto_Image_Pro/DataPro.pl genetumor2.txt >geneTumor.txt
perl /public/home/chshen/auto_Image_Pro/gene_pos_circos_2nd.pl /public/home/chshen/auto_Image_Pro/refGene.txt geneNormal.txt >NormalGene.circos
perl /public/home/chshen/auto_Image_Pro/gene_pos_circos_2nd.pl /public/home/chshen/auto_Image_Pro/refGene.txt geneTumor.txt >TumorGene.circos
perl -lane 'my $a=-$F[3]; print "$F[0]\t$F[1]\t$F[2]\t$a"' NormalGene.circos >NormalGene.circos.p
cp /public/home/chshen/auto_Image_Pro/config_NT_hmgene_g_gte10_2.conf2now ./
circos -conf config_NT_hmgene_g_gte10_2.conf2now
rm -f genenormal.txt genetumor.txt genenormal1.txt genetumor1.txt genenormal2.txt genetumor2.txt geneNormal.txt geneTumor.txt config_NT_hmgene_g_gte10_2.conf2now
