# 1. Install
Download all the files into one folder and run all_in_one.pl using "perl /absolute_path/all_in_one.pl parameter1 parameter2 parameter3 ......".
The users should install the packages used in perl and python programs, such as PerlIO::gzip, Getopt::Long, File::Basename, etc. Also remember to grant executable permissions to some software, or the pipeline will not be able to run. For example:

chmod -R 755 bwa; chmod -R 755 samtools; chmod -R 755 msort

# 2. A Step-to-step protocol of the HIVID2 pipeline 

## 2.1 Step to step tutorial

### stage 1: create the sample list

Manually create a file named total.sample.list(total sample list) should be step1/sample.list. Note that the path in the total.sample.list should be absolute full path and the first four columns should preferably be the same. The header line should be start with #. Below is an example of total.sample.list:

#Sample  FC  Lane  Libray  read_length library_size  fq1  fq2   
SRR12345  SRR12345  SRR12345  SRR12345  110;110 170 /absolute_path/5.fq1.gz /absolute_path/5.fq2.gz  
SRR12346  SRR12346  SRR12346  SRR12346  110;110 170 /absolute_path/6.fq1.gz /absolute_path/6.fq2.gz  
SRR12347  SRR12347  SRR12347  SRR12347  110;110 170 /absolute_path/7.fq1.gz /absolute_path/7.fq2.gz  

### stage 2: run HIVID2 in one single shell script (one-stop pipeline)
perl /absolute_path/all_in_one.pl -o /absolute_path/output_directory -tl /absolute_path/total.sample.list -fa1 /absolute_path/human_ref.fa -fa2 /absolute_path/virus_ref.fa -bin /absolute_path/HIVID2 -c /absolute_path/Config_file

**This program all_in_one.pl is to generate a all-in-one shell script for HIVID2 pipeline**. HIVID2 has 4 steps in stage 2, but using all_in_one.pl, The user can run HIVID2 just by one single shell script.After running all_in_one.pl, each sample will generate a folder in output_directory and there is a shell script named "sampleID_all_in_one.sh" for each sample. Just run it using
"**sh sampleID_all_in_one.sh**". 

Please note that the reference genomes of both human and virus should be indexed by both bwa and soap2 before running sampleID_all_in_one.sh.

**The parameters of all_in_one.pl**
                   
                    -o              <str>           absolute path of output directory
  
                    -tl             <str>           total sample list
  
                    -fa1            <str>           the absolute path of indexed human reference when performing bwa-mem [hg19]
  
                    -fa2            <str>           the absolute path of indexed virus reference when performing bwa-mem [virus]
  
                    -bin            <str>           the absolute path of HIVID2 program
  
                    -c              <str>           the absolute path of configure file for running soap
  
### -c   the Configure file
This configure file difined the indexed referece genomes and alignment parameters used in soap alignment of step3. The users can make their own configure file. But we have involved some configure files which is named as Config* in the same folder of main.pl. Below is the description of the configuration file:  
soap: the path of the soap2 program  
ref_virus: the path of soap2 index of virus reference genome  
ref_human: the path of soap2 index of human reference genome  
insert_sd: the standard deviation of the insert size for the sequencing library  
virus_config: the parameters of soap2 corresponding to different read length; for example, "150;150:-l 50 -v 5 -r 1" means when the read length is 150 bps, then soap2 will use the parameter "-l 50 -v 5 -r 1"; please note that read length is set at sample.list under the folder step1.

The command for creating the configure file is:
python /absolute_path/creat_config.py -soap /absolute_path/soap2 -virus /absolute_path/ref_virus_index -human /absolute_path/ref_human_index -o /absolute_path/Config_file

**soap2 is at soap folder, where 2bwt-builder is used to build the index for the reference genome.** 
  The commmand for running 2bwt-builder is "2bwt-builder ref.fa"
  
## 2.2 Descript of result file and the format

The path of the files of final results:

**The file of final human breakpoint**: step4/*/human/breakpoint/high_confident.*human_bk.final.stp2.uniq2.final

**The file of final virus breakpoint**: step4/*/virus/breakpoint/high_confident.*virus_bk.final.stp2.uniq
 
The low confident breakpoints were stored in files named low_confident.*, please see our paper in published in Bioinformatics for detail.

**Format description of the result file:**

1st column is the number of the chromosome where the breakpoint located.

2nd column is Specific position coordinates

3rd column is the pair amount of left support reads

4th column is the pair amount of right support reads

5th column is the pair amount of discordant support reads

6th column is total number of all support reads

7th column is normalized pair amount of left support reads (normalized_value =support_reads_number / effective_reads_number_of_the_sample)

8th column is normalized pair amount of right support reads (logarithmic) normalized value

9th column is normalized pair amount of discordant support reads

10th column is total number of reads (logarithmic) normalized value

11th column is reads id of left support reads

12th column is reads id of right support reads
  
13th column is reads id of discordant reads supporting the breakpoint

## 2.3 The introduction of main.pl
**main.pl is to generate shell scripts for manualy running 4 steps in stage2 of HIVID2**

Parameters
  
**-o**	   output directory path  
**-l**	   a file containing sample_id, library_id and FC_id  
**-stp**   step number (1/2/3/4)  
**-c**	   parameter configuration file  
**-filter**	   whether to filter the repeated comparison reads. Here, only the repeated comparison reads on the human genome are filtered. The repeated comparison reads on the HBV genome are not filtered. However, in the result, the reads of repeated alignments on the HBV genome will be discarded, and the only aligned reads on the corresponding human genome will be retained.  
**-f**     this parameter is currently useless，please do not use it.

## 2.4 Description of several predefinding files
### (1) -c   the Configure file
This configure file difined the referece genomes and alignment parameters used in step3. The users can make their own configure file. But we have involved some configure files which is named as Config* in the same folder of main.pl. Below is the description of the configuration file:  
soap: the path of the soap2 program  
ref_virus: the path of soap2 index of virus reference genome  
ref_human: the path of soap2 index of human reference genome  
insert_sd: the standard deviation of the insert size for the sequencing library  
virus_config: the parameters of soap2 corresponding to different read length; for example, "150;150:-l 50 -v 5 -r 1" means when the read length is 150 bps, then soap2 will use the parameter "-l 50 -v 5 -r 1"; please note that read length is set at sample.list under the folder step1.

# 3. One demo
A demo has been uploaded. Users can download the file "demo.rar" and unzip it. We have add an file named "used.cml" in each folder. used.cml contains the command lines used in that folder. Please note that users should replace the absolute path of all the files in each script to run the demo. 

# 4. Advanced analysis

After obtaining the integration sites, HIVID2 allows the user to decide whether to automatically perform advanced analysis using the identified virus integrations. 

(1)	Manually seprate result folders of step4 into two groups, For example, tumor and normal, or other user-definednames. If you ran tumor and normal samples in a single run, then you may move each sample (each sample has a folder in step4) into the tumor or normal folder; if you iniatially ran tumor and normal samples seprately during step4, then you can simply use the step4 folder of tumor and normal of each run.

(2)	Run advanced analysis
#First， run Analyse.sh, generatint R scripts and the relevant files.
sh /absolute_path_of_main.pl/advanced_analysis/Analyse.sh /absolute_path/tumor /absolute_path/normal        
#Second, run the generated R scripts
Rscript xxx.R

Note: If you want to get the graph one by one, please separate the script and change parameters. You can also run it line by line, and modify the parameters by yourself. 

# 5. Other tips
(1) In order to help the users to track the data processing, HIVID2 retained some intermediate procedure files during running of the pipeline. It may cause big hard disk consuming when deal with large amount of data such as WGS data. Fortunately, The users can can remove most of intermediate files of previous steps when running step4. When running step4, the user can remove all the files named "*paired.gz" and "*unpaired.gz" in step2, all the files named "*soap.gz" in step2. After completing step4, all the files except the files of final results could be deleted. But before deleting, the users should make sure they don't need them later.

(2) About setting the length in step1/sample.list: the value of read length will not be used in step2 but used in station.sh of step3 for estimating PCR duplication rate. It is OK to set the length based on the raw reads, but it will be better to set the length after running Human_virus_soap.sh in step3 because station.sh of step3 later will directly use this file (Human_*.pair.soap.gz) to estimate the PCR duplication rate of the reads. A simple way is to refer to step3/sample_name/SOAP/Human_*.pair.soap.gz for estimating read length. Length of each read is given in the soap file. You can use a length with the highest frequency as the length to set in sample.list as the length of each read will be not the same. For example, if there are five reads in total and the read length after trimmomatics are: 100, 95, 100, 100, 100, then 100 should be set as the length in sample.list.
    

Alternatively, users can set the read length in sample.list after completing step2 (trimmomatics, quality control). Please note that the length of each read will be not the same after trimmomtics. So You can use a length with highest frequency as the length to set in sample.list. For example, if there are five reads in total and the read length after trimmomatics are: 100, 95, 100, 100, 100, then 100 should be set as the length in sample.list. You can refer to the file "*.trimmo.paired.gz" (fq1 or fq2) in step2 for estimating the reads length. 
 
(3) There is a file named "tfbsConsSites.txt" in the advanced analysis. This file cannot be uploaded onto github due to the size limitation. But the user could download this file from Table browser of UCSC.

(4) HIVID2 works quite well for virus-capture sequencing data. For WGS data, sometimes the used memory might be too large. In this case, the users may need to separate the fastq data into several parts before input into HIVID2 for step1,step2 and step3; then the users can merge the data of step3 for the separated parts to run step4. For WGS data, the users could alternatively first remove human reads or HBV reads before running HIVID2. 

(5) The programs of soap2 are also uploaded in this repository. The method of buliding index of soap2 is as following:
    2bwt-builder xx.fa
 
(6) It should be noted that there are a file named "ref.list" in the same folder of main.pl. "ref.list" must contain all the ID of reference genomes used in the sequence alignment of step3 and step4, or the user will get error or uncompleted results in *human_bk.final.stp2.uniq2.final during the procedure of deep removing PCR-duplications in step4. We have involved some predefined reference names in ref.list, but the users should add the references names used in their own experiments. In the ref.list, each ID should be followed by an underline, for example "chr1_".

(7) The workflow of HIVID2 includes four sub-steps, if you want to manually run the 4 steps one by one, please refer to the file named "HIVID2_manunally_run_four_steps.docx".

(8) For whether to use normalized number of support reads: as the normalization value might not be very precise due to the difficulty of estimation the PCR duplication for raw reads, the user could choose to use absolute number of support reads or normalized number of support reads depending  on the real situation.

# 6. Citation
Xi Zeng, Linghao Zhao, Chenhang Shen, Yi Zhou, Guoliang Li, Wing-Kin Sung, HIVID2: an accurate tool to detect virus integrations in the host genome, Bioinformatics, Volume 37, Issue 13, 1 July 2021, Pages 1821–1827, https://doi.org/10.1093/bioinformatics/btab031
