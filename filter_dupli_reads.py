# -*- coding: utf-8 -*-
#！/usr/bin/python
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-se", type=str)
parser.add_argument("-R1", type=str)
parser.add_argument("-R2", type=str)
parser.add_argument("-out1", type=str)
parser.add_argument("-out2", type=str)
args = parser.parse_args()
disc_read=open(args.se,'r')
#disc_read=open('/public/home/wyyy/HPV/HIVID2.1/step3_update/20C257574/station_pair_end/20C257574_se_se','r')
disc_readID=[]
for line in disc_read.readlines():
    id=line.split('\t')[0]
    disc_readID.append(id)
disc_readID=set(disc_readID)
n=0
allread1=[]
allread2=[]
with open(args.R1,'r') as handle:
#with open('/public/home/wyyy/HPV/HIVID2.1/step2/20C257574/20C257574_R302_CapNGS.clean_R1.fq.trimmo.paired','r') as handle:
    fq=SeqIO.parse(handle,'fastq')
    for reads in fq:
        for i in disc_readID:
            if reads.id==i.split('/')[0] and i.split('/')[1]=='1':
                allread1.append(reads)
#SeqIO.write(allread1,'/public/home/wyyy/HPV/HIVID2.1/step3_update/20C257574/20C257574.disc.R1.fastq','fastq')
SeqIO.write(allread1,args.out1,'fastq')
with open(args.R2,'r') as handle:
#with open('/public/home/wyyy/HPV/HIVID2.1/step2/20C257574/20C257574_R302_CapNGS.clean_R2.fq.trimmo.paired','r') as handle:
    fq=SeqIO.parse(handle,'fastq')
    for reads in fq:
        for i in disc_readID:
            if reads.id==i.split('/')[0] and i.split('/')[1]=='2':
                allread2.append(reads)
SeqIO.write(allread1,args.out2,'fastq')
#SeqIO.write(allread2,'/public/home/wyyy/HPV/HIVID2.1/step3_update/20C257574/20C257574.disc.R2.fastq','fastq')
