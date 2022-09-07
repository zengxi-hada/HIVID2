mkdir $1/step3/$2/station_pair_end_new
mkdir $1/step3/$2/reads_assemble_pair-end_new
perl $3/new_HBV_human_soap.pl -reads_assembly $3/overlap_pair_trim.new -margefa_qua $3/margefa_qua.pl -hp $1/step3/$2/SOAP/Human_$2.pair.soap.gz -hs $1/step3/$2/SOAP/Human_$2.single.soap.gz -hu $1/step3/$2/SOAP/Human_$2.unmap.soap.gz -bp $1/step3/$2/SOAP/virus_$2.pair.soap.gz -bs $1/step3/$2/SOAP/virus_$2.single.soap.gz -bu $1/step3/$2/SOAP/virus_$2.unmap.soap.gz -pe $1/step3/$2/station_pair_end_new/$2_pe_pe.gz -se $1/step3/$2/station_pair_end_new/$2_se_se.gz -sb $1/step3/$2/station_pair_end_new/$2_virus_un.gz -sh $1/step3/$2/station_pair_end_new/$2_Human_un.gz -un $1/step3/$2/station_pair_end_new/$2_un_un.gz  -stat $1/step3/$2/$2.stat -f1 $1/step2/$2/$filename1.trimmo.paired.gz -f2 $1/step2/$2/$filename2.trimmo.paired.gz -o $1/step3/$2/reads_assemble_pair-end_new
gunzip $1/step3/$2/station_pair_end_new/$2_se_se.gz
filename=`sed -n "2,2p" $1/step1/sample.list`
filename1=`echo $filename | awk '{print $7}'`;filename1=${filename1##*/};filename1=${filename1%.*}
filename2=`echo $filename | awk '{print $8}'`;filename2=${filename2##*/};filename2=${filename2%.*}
gunzip $1/step2/$2/$filename1.trimmo.paired.gz
gunzip $1/step2/$2/$filename2.trimmo.paired.gz
gunzip $1/step3/$2/SOAP/Human_$2.single.soap.gz
gunzip $1/step3/$2/SOAP/virus_$2.single.soap.gz
python $3/filter_dupli_reads.py -se $1/step3/$2/station_pair_end_new/$2_se_se -R1 $1/step2/$2/$filename1.trimmo.paired -R2 $1/step2/$2/$filename2.trimmo.paired -out1 $1/step4/$2/fq/$2.disc.R1.fastq -out2 $1/step4/$2/fq/$2.disc.R2.fastq
python $3/stat_reads_length.py -R1 $1/step4/$2/fq/$2.disc.R1.fastq -R2 $1/step4/$2/fq/$2.disc.R2.fastq -out1 $1/step4/$2/fq/$2_R1.txt -out2 $1/step4/$2/fq/$2_R2.txt
gzip $1/step4/$2/fq/$2.disc.R1.fastq
gzip $1/step4/$2/fq/$2.disc.R2.fastq
bwa mem $4 $1/step4/$2/fq/$2.disc.R1.fastq.gz $1/step4/$2/fq/$2.disc.R2.fastq.gz > $1/step4/$2/human/human_$2.disc.sam
samtools view -h -q 9 $1/step4/$2/human/human_$2.disc.sam > $1/step4/$2/human/human_$2.disc.filter.sam
bwa mem $5 $1/step4/$2/fq/$2.disc.R1.fastq.gz $1/step4/$2/fq/$2.disc.R2.fastq.gz > $1/step4/$2/virus/virus_$2.disc.sam
samtools view -h -q 9 $1/step4/$2/virus/virus_$2.disc.sam > $1/step4/$2/virus/virus_$2.disc.filter.sam
python $3/se_se.v2.py -se $1/step3/$2/station_pair_end_new/$2_se_se -humansoap $1/step3/$2/SOAP/Human_$2.single.soap -virussoap $1/step3/$2/SOAP/virus_$2.single.soap -r1 $1/step4/$2/fq/$2_R1.txt -r2 $1/step4/$2/fq/$2_R2.txt -o $1/step3/$2/station_pair_end_new/$2_se_se_new
python $3/filter_se_se_reads.py -se $1/step3/$2/station_pair_end_new/$2_se_se_new -humansam $1/step4/$2/human/human_$2.disc.filter.sam -virussam $1/step4/$2/virus/virus_$2.disc.filter.sam -o $1/step3/$2/station_pair_end_new/$2_se_se_new_update
python $3/filter_se_se_reads_HBV.py -se $1/step3/$2/station_pair_end_new/$2_se_se_new -humansam $1/step4/$2/human/human_$2.disc.filter.sam -virussam $1/step4/$2/virus/virus_$2.disc.filter.sam -o $1/step3/$2/station_pair_end_new/$2_se_se_new_update_forHBV
gzip $1/step3/$2/station_pair_end_new/$2_se_se_new
gzip $1/step3/$2/station_pair_end_new/$2_se_se_new_update
gzip $1/step3/$2/station_pair_end_new/$2_se_se_new_update_forHBV
perl $3/breakpiont_discordant-rd_v2.pl -disc_read $1/step3/$2/station_pair_end_new/$2_se_se_new_update.gz -bk $1/step4/$2/human/breakpoint/$2_R302_CapNGS_human_bk.final.stp2.uniq -len 250 -ods $1/step4/$2/human/breakpoint/human_bk.discordant_read_new.txt1.2 -obk $1/step4/$2/human/breakpoint/$2_R302_CapNGS_human_bk.final.stp2.uniq2
perl $3/breakpiont_discordant_forHBV-rd_v2.pl -disc_read $1/step3/$2/station_pair_end_new/$2_se_se_new_update_forHBV.gz -bk $1/step4/$2/virus/breakpoint/high_confident_$2_R302_CapNGS_virus_bk.final.stp2.uniq -len 50 -ods $1/step4/$2/virus/breakpoint/virus_bk.discordant_read_new.txt1.2 -obk $1/step4/$2/virus/breakpoint/high_confident_$2_R302_CapNGS_virus_bk.final.stp2.uniq2
sed -ir "s#nan#0#g" $1/step4/$2/human/breakpoint/$2_R302_CapNGS_human_bk.final.stp2.uniq2
sed -ir "s#nan#0#g" $1/step4/$2/virus/breakpoint/high_confident_$2_R302_CapNGS_virus_bk.final.stp2.uniq2
python $3/try-fast.v3.py3.py -fq1 $1/step2/$2/$2_R302_CapNGS.clean_R1.fq.trimmo.paired1.gz -fq2 $1/step2/$2/$2_R302_CapNGS.clean_R2.fq.trimmo.paired1.gz -i $1/step4/$2/human/breakpoint/$2_R302_CapNGS_human_bk.final.stp2.uniq2 -o $1/step4/$2/human/breakpoint/high_confident_$2_R302_CapNGS_human_bk.final.stp2.uniq2.final -id $2 -ref $3/ref.list
gzip $1/step3/$2/station_pair_end_new/$2_se_se
gzip $1/step2/$2/$2_R302_CapNGS.clean_R1.fq.trimmo.paired
gzip $1/step2/$2/$2_R302_CapNGS.clean_R2.fq.trimmo.paired
gzip $1/step3/$2/SOAP/Human_$2.single.soap
gzip $1/step3/$2/SOAP/virus_$2.single.soap
