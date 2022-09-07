# -*- coding: utf-8 -*-
import gzip
import os
from Bio import SeqIO
from tempfile import TemporaryFile
import argparse
import sys
import re
reload(sys) 
sys.setdefaultencoding('utf-8')
sys.setrecursionlimit(10000)

## author
#Zhou Yi

bin_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
temp = TemporaryFile()

cmd_parser = argparse.ArgumentParser(description='re-evaluate uniq_read in support_read')
cmd_parser.add_argument('-fq1', help='fastq1 after remove adapter')
cmd_parser.add_argument('-fq2', help='fastq2 after remove adapter')
cmd_parser.add_argument('-i', help='human-uniq infile')
cmd_parser.add_argument('-o', help='outfile, directory and name')
cmd_parser.add_argument('-erate',type=float,default=0.04,help='error rate of sequencing for allowing mismatch')
cmd_parser.add_argument('-ref',default="./ref.list",help='ref genome ID list')
cmd_parser.add_argument('-id',help='sample name')
cmd_args = cmd_parser.parse_args()

#互补序列  不需要反向
def revseq(seq):
	seq=seq.replace("A","t").replace("T","a").replace("G","c").replace("C","g")
	seq=seq.upper()  #[::-1]
	return seq

#两个序列比较
def compared_readseq(seq1,seq2,mismatch_percent=0.04):
	if len(seq1) != len(seq2):
		return 0  #非重复
	else:
		mismatch=mismatch_percent*len(seq2)
		count_mis=0
		for i in range(0,len(seq1)):
			if seq1[i] != seq2[i]:
				count_mis+=1
		if count_mis <= (mismatch):
			return 1  #重复的
		else:
			return 0  #非重复

#判断seq in seq_list
def compared_readseq_list(seq,list1):
	for seq_tmp in list1:
		if compared_readseq(seq,seq_tmp):
			return 1  #1是重复的
			break
		else:
			continue
	return 0

#获得fq文件里某一id的read的seq,返回字典
def get_seq(file1,file2,id_list):
	dic1={}
	dic2={}
	tmp=cmd_args.id+".list"
	tmp_o1=cmd_args.id+".tmp1.fq"
	tmp_o2=cmd_args.id+".tmp2.fq"
	with open(tmp,"w") as out1:
		out1.write("\n".join(id_list)+"\n")
	cmd1="zcat " +file1+ " | " +bin_dir+ "/seqkit grep -f " +tmp+ " -o " +tmp_o1
	cmd2="zcat " +file2+ " | " +bin_dir+ "/seqkit grep -f " +tmp+ " -o " +tmp_o2
	os.system(cmd1)
	os.system(cmd2)
	os.remove(tmp)
	#handle=os.system("zcat file1 | ~/tools/seqkit grep -f temp")
	#handle = gzip.open("tmp.fq", "r")
	#record = (r1 for r1 in SeqIO.parse(handle, "fastq") if r1.id in id_list)
	for r1 in SeqIO.parse(tmp_o1, "fastq"):
		dic1[str(r1.id)]=str(r1.seq)
		#print str(r1.id)
	for r2 in SeqIO.parse(tmp_o2, "fastq"):
		dic2[str(r2.id)]=str(r2.seq)
		#print str(r2.id)
	os.remove(tmp_o1)
	os.remove(tmp_o2)
	return dic1,dic2

#从dic获取uniq list id的列表
def get_uniq(dic1,dic2):
	id_list=dic1.keys()
	if 'uniq_list' not in globals():
		global uniq_list
		uniq_list=[]
	if len(id_list)==0:
#		print "test1"
		return []
	elif len(id_list)==1:
#		print "test2"
		uniq_list.append(id_list[0])
		#print uniq_list
		return uniq_list
	else:
#		print "test3"
		seq_list1=dic1[id_list[0]]
		seq_list2=dic2[id_list[0]]
		dic1.pop(id_list[0])
		dic2.pop(id_list[0])
		seq_list1_tmp=dic1.values()
		seq_list2_tmp=dic2.values()

		a1=compared_readseq_list(seq_list1,seq_list1_tmp)
		b1=compared_readseq_list(seq_list2,seq_list2_tmp)

		a2=compared_readseq_list(seq_list1,seq_list2_tmp)
		b2=compared_readseq_list(seq_list2,seq_list1_tmp)

		a3=compared_readseq_list(revseq(seq_list1),seq_list1_tmp)
		b3=compared_readseq_list(revseq(seq_list2),seq_list2_tmp)

		a4=compared_readseq_list(revseq(seq_list1),seq_list2_tmp)
		b4=compared_readseq_list(revseq(seq_list2),seq_list1_tmp)
		if ((a1==1 and b1 ==1) or (a2==1 and b2 ==1) or (a3==1 and b3 ==1) or (a4==1 and b4 ==1))==0:
			uniq_list.append(id_list[0])
			#print uniq_list
		return get_uniq(dic1,dic2)


#对一个read,提取原始的readID
def clean_readID(read_id,ref_list):
	if read_id[-2:]=="_1" or read_id[-2:]=="_2":
		read_id=read_id[0:-2]
	elif read_id[-2:]=="/1" or read_id[-2:]=="/2":
		read_id=read_id[0:-2]
	else:
		read_id=read_id
	read_id=read_id.replace("left_","").replace("trim_pe#","").replace("trim_se#","").replace("unmap_unmap_","").replace("unmap_unmap","").replace("unmap_","")
	for i in ref_list:
		read_id=read_id.replace(i,"")
#	read_id=read_id.replace("_","")
		
	return read_id


#main
ref_list=[]
with open (cmd_args.ref,"r") as reflist:
	for line in reflist:
		line=line.strip("\n")
		line=re.sub(r'\s+', '', line)
		ref_list.append(line)

with open (cmd_args.i,"r") as infile:
	with open (cmd_args.o,"w") as outfile:
		for line in infile:
			line=line.strip("\n")
##			print line
			line=line.strip("\n").split("\t")
			if line[0] !="ref":
				#left_reads_ID
				raw_id_list1={}
				raw_id_list2={}
				raw_id_list3={}
				if line[-3] == "0":
					uniq_read1=0
					uniq_id_list1=["0"]
					raw_id_list1["0"]="0"
				else:
					left_reads_ID=line[-3].split(",")
					if len(left_reads_ID) ==0:
						print "error,file without left_reads_ID.\n"
						break
					elif len(left_reads_ID) ==1:
						uniq_id_list1=[left_reads_ID[0]]
						uniq_read1=1
						raw_id_list1[left_reads_ID[0]]=left_reads_ID[0]
					else:
						id_list1=[]
						for ids1 in left_reads_ID:
							clean_id1=clean_readID(ids1,ref_list)
							id_list1.append(clean_id1)
							raw_id_list1[clean_id1]=ids1
						#print id_list1
						uniq_list=[]
						dic1,dic2=get_seq(cmd_args.fq1,cmd_args.fq2,id_list1)
						uniq_id_list1=get_uniq(dic1,dic2)
						uniq_read1=len(uniq_id_list1)
				#right_reads_ID
				if line[-2] == "0":
					uniq_read2=0
					uniq_id_list2=["0"]
					raw_id_list2["0"]="0"
				else:
					right_reads_ID=line[-2].split(",")
					if len(right_reads_ID) == 0:
						print "error,file without right_reads_ID.\n"
						break
					elif len(right_reads_ID) == 1:
						uniq_id_list2=[right_reads_ID[0]]
						uniq_read2=1
						raw_id_list2[right_reads_ID[0]]=right_reads_ID[0]
					else:
						id_list2=[]
						for ids2 in right_reads_ID:
							clean_id2=clean_readID(ids2,ref_list)
							id_list2.append(clean_id2)
							raw_id_list2[clean_id2]=ids2
						#print id_list2
						uniq_list=[]
						dic1,dic2=get_seq(cmd_args.fq1,cmd_args.fq2,id_list2)
						uniq_id_list2=get_uniq(dic1,dic2)
						uniq_read2=len(uniq_id_list2)
				#discordant_reads_ID (zengxi)
				if line[-1] == "0":
					uniq_read3=0
					uniq_id_list3=["0"]
					raw_id_list3["0"]="0"
##					print "zero";
				else:
					disc_reads_ID=line[-1].split(",")					#discrodant reads ID
					if len(disc_reads_ID) == 0:
						print "error,file without right_reads_ID.\n"
						break
##						print disc_reads_ID
					elif len(disc_reads_ID) == 1:
						uniq_id_list3=[disc_reads_ID[0]]
						uniq_read3=1
						raw_id_list3[disc_reads_ID[0]]=disc_reads_ID[0]
##						print disc_reads_ID[0]
					else:
						id_list3=[]
						for ids3 in disc_reads_ID:
							clean_id3=clean_readID(ids3,ref_list)
							id_list3.append(clean_id3)
							raw_id_list3[clean_id3]=ids3
##							print ids3
						#print id_list3
						uniq_list=[]
						dic1,dic2=get_seq(cmd_args.fq1,cmd_args.fq2,id_list3)
						uniq_id_list3=get_uniq(dic1,dic2)
						uniq_read3=len(uniq_id_list3)
	
				#print uniq_id_list1
				output1="\t".join(line[0:2])
				left_support=str(uniq_read1)
				right_support=str(uniq_read2)
				disc_support=str(uniq_read3)
				total_support=str(uniq_read1+uniq_read2+uniq_read3)
				norm_left=str(round(uniq_read1/float(line[2])*float(line[6]),3)) if uniq_read1>0 else "0"
				norm_right=str(round(uniq_read2/float(line[3])*float(line[7]),3)) if uniq_read2>0 else "0"
				norm_disc=str(round(uniq_read3/float(line[4])*float(line[8]),3)) if uniq_read3>0 else "0"
				norm_sum=str(round(float(norm_left)+float(norm_right)+float(norm_disc),3))
				left_reads_ID=",".join(raw_id_list1[a] for a in uniq_id_list1)
				right_reads_ID=",".join(raw_id_list2[b] for b in uniq_id_list2)
				disc_reads_ID=",".join(raw_id_list3[c] for c in uniq_id_list3)
				#outfile.write("\t".join(line[:-2])+"\t"+str(uniq_read1)+"\t"+str(uniq_read2)+"\t"+str(uniq_read1+uniq_read2)+"\t"+",".join(uniq_id_list1)+"\t"+",".join(uniq_id_list2)+"\n")
				outfile.write(output1+"\t"+left_support+"\t"+right_support+"\t"+disc_support+"\t"+total_support+"\t"+norm_left+"\t"+norm_right+"\t"+norm_disc+"\t"+norm_sum+"\t"+left_reads_ID+"\t"+right_reads_ID+"\t"+disc_reads_ID+"\n")
			else:
				outfile.write("\t".join(line)+"\n")
				#outfile.write("\t".join(line[:-2])+"\tleft_uniq\tright_uniq\tTotal_uniq\tleft_reade_new\tright_reads_new\n")
