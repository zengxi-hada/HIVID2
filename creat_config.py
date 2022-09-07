# -*- coding: utf-8 -*-
#！/usr/bin/python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-soap", type=str)
parser.add_argument("-virus", type=str)
parser.add_argument("-human", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()
output=open(args.o,'a')
print('>soap=',args.soap,'\n','>ref_virus=',args.virus,'\n','>ref_human=',args.human,file=output,sep='')
print('>insert_sd=-30/+30\n>virus_config=\n29;29:-v 1 -r 1\n44;44:-v 2 -r 1\n45;45:-v 2 -r 1\n35;44:-v 2 -r 1\n51;51:-v 2 -r 1\n50;50:-v 2 -r 1\n49;49:-v 2 -r 1\n60;60:-v 2 -r 1\n75;75:-l 40 -v 4 -r 1\n76;76:-l 40 -v 4 -r 1\n73;75:-l 40 -v 4 -r 1\n72;75:-l 40 -v 4 -r 1\n90;90:-l 40 -v 5 -r 1\n100;100:-l 40 -v 5 -r 1\n141;141:-l 50 -v 5 -r 1\n150;150:-l 50 -v 5 -r 1\n>Human_config=\n29;29:-v 1 -r 1\n44;44:-v 2 -r 1\n45;45:-v 2 -r 1\n35;44:-v 2 -r 1\n51;51:-v 2 -r 1\n50;50:-v 2 -r 1\n49;49:-v 2 -r 1\n60;60:-v 2 -r 1\n75;75:-l 40 -v 4 -r 1\n76;76:-l 40 -v 4 -r 1\n73;75:-l 40 -v 4 -r 1\n72;75:-l 40 -v 4 -r 1\n90;90:-l 40 -v 5 -r 1\n100;100:-l 40 -v 5 -r 1\n141;141:-l 50 -v 5 -r 1\n150;150:-l 50 -v 5 -r 1',file=output)
