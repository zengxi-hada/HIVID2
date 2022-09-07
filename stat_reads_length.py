# -*- coding: utf-8 -*-
#！/usr/bin/python
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser()
parser.add_argument("-R1", type=str)
parser.add_argument("-R2", type=str)
parser.add_argument("-out1", type=str)
parser.add_argument("-out2", type=str)
args = parser.parse_args()
output1=open(args.out1,'a')
output2=open(args.out2,'a')
with open(args.R1,'r') as handle:
    fq=SeqIO.parse(handle,'fastq')
    for reads in fq:
        print(reads.id,len(reads.seq),sep='\t',file=output1)
with open(args.R2,'r') as handle:
    fq=SeqIO.parse(handle,'fastq')
    for reads in fq:
        print(reads.id,len(reads.seq),sep='\t',file=output2)
