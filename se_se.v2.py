# -*- coding: utf-8 -*-
#！/usr/bin/python
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-se", type=str)
parser.add_argument("-humansoap", type=str)
parser.add_argument("-virussoap", type=str)
parser.add_argument("-r1", type=str)
parser.add_argument("-r2", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()
i=0
file = gzip.open(args.se,'rt')
output=open(args.o,'a')
printline=[]
human_align={}
virus_align={}

with open(args.humansoap,'r') as strandfile1:
    for reads in strandfile1:
        readsid=reads.split('\t')[0]
        readschr=reads.split('\t')[7]
        readspos=reads.split('\t')[8]
        key=readsid+','+readschr+','+readspos
        human_align[key]=reads.split('\t')[6]
strandfile1.close()
with open(args.virussoap,'r') as strandfile2:
    for reads in strandfile2:
        readsid=reads.split('\t')[0]
        readschr=reads.split('\t')[7]
        readspos=reads.split('\t')[8]
        key=readsid+','+readschr+','+readspos
        virus_align[key]=reads.split('\t')[6]
strandfile2.close()
for line in file:
    id=line.split('\t')[0].split('/')[0]
    read=line.split('\t')[0].split('/')[1]
    chr=line.split('\t')[1]
    pos=line.strip().split('\t')[2]
    key=id+','+chr+','+pos
    if read=='1':
        with open(args.r1,'r') as fq1:
            for reads in fq1:
                readsid=reads.split('\t')[0]
                if id==readsid:
                    read1=line.strip()
                    len1=reads.strip().split('\t')[1]
        fq1.close()
        if line.split('\t')[1][:3]=='chr':
            strand1=human_align[key]
            line1=read1+'\t'+strand1
            printline.append(line1)
        else:
            strand1=virus_align[key]
            line1=read1+'\t'+strand1
            printline.append(line1)
        i=i+1
    if read=='2':
        with open(args.r2,'r') as fq2:
            for reads in fq2:
                readsid=reads.split('\t')[0]
                if id==readsid:
                    read2=line.strip()
                    len2=reads.strip().split('\t')[1]
        fq2.close()
        if line.split('\t')[1][:3]=='chr':
            strand2=human_align[key]
            line2=read2+'\t'+strand2
            printline.append(line2)
        else:
            strand2=virus_align[key]
            line2=read2+'\t'+strand2
            printline.append(line2)
        i=i+1
    if i==2:
        l1=printline[0]+'\t'+len1+'\t'+len2
        l2=printline[1]+'\t'+len1+'\t'+len2
        print(l1,file=output)
        print(l2,file=output)
        i=0
        strand1=''
        len1=''
        strand2=''
        len2=''
        printline=[]
