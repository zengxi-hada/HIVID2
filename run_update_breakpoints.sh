filename=`sed -n "2,2p" $1/step1/sample.list`
#filename1=`echo $filename | awk '{print $7}'`;filename1=${filename1##*/};filename1=${filename1%.*}
#filename2=`echo $filename | awk '{print $8}'`;filename2=${filename2##*/};filename2=${filename2%.*}
filename1=`echo $filename | awk '{print $7}'`;filename1=${filename1##*/};filename1=${filename1//.gz/}
filename2=`echo $filename | awk '{print $8}'`;filename2=${filename2##*/};filename2=${filename2//.gz/}

mkdir $1/step3/$2/station_pair_end_new
mkdir $1/step3/$2/reads_assemble_pair-end_new
perl $3/new_HBV_human_soap.pl -reads_assembly $3/overlap_pair_trim.new -margefa_qua $3/margefa_qua.pl -hp $1/step3/$2/SOAP/Human_$2.pair.soap -hs $1/step3/$2/SOAP/Human_$2.single.soap -hu $1/step3/$2/SOAP/Human_$2.unmap.soap -bp $1/step3/$2/SOAP/virus_$2.pair.soap -bs $1/step3/$2/SOAP/virus_$2.single.soap -bu $1/step3/$2/SOAP/virus_$2.unmap.soap -pe $1/step3/$2/station_pair_end_new/$2_pe_pe.gz -se $1/step3/$2/station_pair_end_new/$2_se_se.gz -sb $1/step3/$2/station_pair_end_new/$2_virus_un.gz -sh $1/step3/$2/station_pair_end_new/$2_Human_un.gz -un $1/step3/$2/station_pair_end_new/$2_un_un.gz  -stat $1/step3/$2/$2.stat -f1 $1/step2/$2/$filename1.trimmo.paired.gz -f2 $1/step2/$2/$filename2.trimmo.paired.gz -o $1/step3/$2/reads_assemble_pair-end_new

python $3/filter_dupli_reads.py -se $1/step3/$2/station_pair_end_new/$2_se_se.gz -R1 $1/step2/$2/$filename1.trimmo.paired.gz -R2 $1/step2/$2/$filename2.trimmo.paired.gz -out1 $1/step4/$2/fq/$2.disc.R1.fastq -out2 $1/step4/$2/fq/$2.disc.R2.fastq
python $3/stat_reads_length.py -R1 $1/step4/$2/fq/$2.disc.R1.fastq -R2 $1/step4/$2/fq/$2.disc.R2.fastq -out1 $1/step4/$2/fq/$2_R1.txt -out2 $1/step4/$2/fq/$2_R2.txt
gzip -f $1/step4/$2/fq/$2.disc.R1.fastq
gzip -f $1/step4/$2/fq/$2.disc.R2.fastq
bwa mem $4 $1/step4/$2/fq/$2.disc.R1.fastq.gz $1/step4/$2/fq/$2.disc.R2.fastq.gz > $1/step4/$2/human/human_$2.disc.sam
samtools view -h -q 9 $1/step4/$2/human/human_$2.disc.sam > $1/step4/$2/human/human_$2.disc.filter.sam
bwa mem $5 $1/step4/$2/fq/$2.disc.R1.fastq.gz $1/step4/$2/fq/$2.disc.R2.fastq.gz > $1/step4/$2/virus/virus_$2.disc.sam
samtools view -h -q 9 $1/step4/$2/virus/virus_$2.disc.sam > $1/step4/$2/virus/virus_$2.disc.filter.sam
python $3/se_se.v2.py -se $1/step3/$2/station_pair_end_new/$2_se_se.gz -humansoap $1/step3/$2/SOAP/Human_$2.single.soap -virussoap $1/step3/$2/SOAP/virus_$2.single.soap -r1 $1/step4/$2/fq/$2_R1.txt -r2 $1/step4/$2/fq/$2_R2.txt -o $1/step3/$2/station_pair_end_new/$2_se_se_new
python $3/filter_se_se_reads.py -se $1/step3/$2/station_pair_end_new/$2_se_se_new -humansam $1/step4/$2/human/human_$2.disc.filter.sam -virussam $1/step4/$2/virus/virus_$2.disc.filter.sam -o $1/step3/$2/station_pair_end_new/$2_se_se_new_update
python $3/filter_se_se_reads_HBV.py -se $1/step3/$2/station_pair_end_new/$2_se_se_new -humansam $1/step4/$2/human/human_$2.disc.filter.sam -virussam $1/step4/$2/virus/virus_$2.disc.filter.sam -o $1/step3/$2/station_pair_end_new/$2_se_se_new_update_forHBV
#rm -f $1/step3/$2/station_pair_end_new/$2_se_se_new
#rm -f $1/step3/$2/station_pair_end_new/$2_se_se_new_update
#rm -f $1/step3/$2/station_pair_end_new/$2_se_se_new_update_forHBV
perl -f $3/breakpiont_discordant-rd_v2.pl -disc_read $1/step3/$2/station_pair_end_new/$2_se_se_new_update -bk $1/step4/$2/human/breakpoint/$2_human_bk.final.stp2.uniq -len 250 -ods $1/step4/$2/human/breakpoint/human_bk.discordant_read_new.txt1.2 -obk $1/step4/$2/human/breakpoint/$2_human_bk.final.stp2.uniq2
perl $3/breakpiont_discordant_forHBV-rd_v2.pl -disc_read $1/step3/$2/station_pair_end_new/$2_se_se_new_update_forHBV -bk $1/step4/$2/virus/breakpoint/high_confident_$2_virus_bk.final.stp2.uniq -len 50 -ods $1/step4/$2/virus/breakpoint/virus_bk.discordant_read_new.txt1.2 -obk $1/step4/$2/virus/breakpoint/high_confident_$2_virus_bk.final.stp2.uniq2
sed -ir "s#nan#0#g" $1/step4/$2/human/breakpoint/$2_human_bk.final.stp2.uniq2
sed -ir "s#nan#0#g" $1/step4/$2/virus/breakpoint/high_confident_$2_virus_bk.final.stp2.uniq2
python $3/try-fast.v3.py3.py -fq1 $1/step2/$2/$filename1.trimmo.paired1.gz -fq2 $1/step2/$2/$filename2.trimmo.paired1.gz -i $1/step4/$2/human/breakpoint/$2_human_bk.final.stp2.uniq2 -o $1/step4/$2/human/breakpoint/high_confident_$2_human_bk.final.stp2.uniq2.final -id $2 -ref $3/ref.list

rm -f $1/step2/$2/*fq*paired*
rm -f $1/step3/$2/SOAP/*soap
rm -f $1/step3/$2/station_pair_end_new/$2_se_se
rm -f $1/step3/$2/station_pair_end_new/$2_se_se_new
rm -f $1/step3/$2/station_pair_end_new/$2_se_se_new_update
rm -f $1/step3/$2/station_pair_end_new/$2_se_se_new_update_forHBV
rm -f $1/step3/$2/reads_assemble_pair-end_new/*/*paired
rm -f $1/step3/$2/reads_assemble_pair-end_new/*/*fa
rm -f $1/step3/$2/reads_assemble_single-end_new/*/*paired
rm -f $1/step4/$2/*/*bam
rm -f $1/step4/$2/*/*sam
rm -f $1/step4/$2/fq/*fq.gz
rm -f $1/step4/$2/fq/*fastq.gz
rm -f $1/step4/$2/fq/*fq

echo \"All work has done!\" >&2\n
