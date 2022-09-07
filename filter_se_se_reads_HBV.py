# -*- coding: utf-8 -*-
#！/usr/bin/python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-se", type=str)
parser.add_argument("-humansam", type=str)
parser.add_argument("-virussam", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()
se=open(args.se,'r')
human=open(args.humansam,'r')
virus=open(args.virussam,'r')
output=open(args.o,'a')
hi=0
vi=0
line=0
for read in se.readlines():
    id=read.split('\t')[0].split('/')[0]
    line=line+1
    if read.split('\t')[1][0]=='c':
        for align in human:
            if align[0]!='@' and id==align.split('\t')[0]:
                humanread=read.strip()
                hi=hi+1
    else:
        for align in virus:
            if align[0]!='@' and id==align.split('\t')[0]:
                virusread=read.strip()
                vi=vi+1
    if line==2:
        if hi==vi==1:
            print(virusread,file=output)
            print(humanread,file=output)
            line=0
            hi=0
            vi=0
        elif hi>1 and vi==1:
            print(humanread)
            print(virusread)
            line=0
            hi=0
            vi=0
        elif hi==1 and vi>1:
            print(humanread)
            print(virusread)
            line=0 
            hi=0
            vi=0
        else:
            line=0
            hi=0
            vi=0
