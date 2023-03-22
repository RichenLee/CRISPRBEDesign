#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import regex as re
import logging
import argparse
import datetime
import sys
import platform


__version__ = "1.0.0"

logging.basicConfig(
					 format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
					 datefmt='%a, %d %b %Y %H:%M:%S',
					 stream=sys.stderr,
					 filemode="w"
					 )
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
error   = logger.critical
warn    = logger.warning
debug   = logger.debug
info    = logger.info

base_dict={'R':'AG','Y':'CT','M':'AC','K':'GT','S':'GC','W':'AT',
			'H':'ATC','B':'GTC','V':'GAC','D':'GAT','N':'ATCG',
			'A':'A','T':'T','C':'C','G':'G'}

def tell_plat():
	if platform.system().lower() == 'windows':
		flag=True
	elif platform.system().lower() == 'linux':
		flag=False
	return flag

def get_header():
	return (r'''
				CRISPRBEDesign                                
 ---a program for base editing sgRNA design                                                   
 Version V1.0.0                                                    
 LAST REVISED: 2023.03.20                                          
''')

def find_all(string, sub):
	start = 0
	pos = []
	while True:
		start = string.find(sub, start)
		if start == -1:
			return pos
		pos.append(start)
		start += len(sub)
	return pos

def tell_window(sequence,base,window_size,window_start):
	flag=False
	loc=''
	window_start=window_start-1
	window_sequence=sequence[window_start:(window_start+window_size)]
	if base in window_sequence:
		flag=True
	if flag:
		loc_list=find_all(window_sequence,base)
		loc=','.join([str(i+window_start+1) for i in loc_list])
	return flag,loc

def count_GC(sequence):
	G=sequence.count('G')
	C=sequence.count('C')
	GC=(G+C)/len(sequence)
	return str(round(GC,2))

def Fasta_reverse(sequence):
	#将序列进行方向互补
	sequence=sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	sequence = sequence.upper()
	return sequence[::-1]

def mode(PAM_end,PAM,spacer_length):
	word=''
	reverse_word=''
	for base in PAM:
		word+='[%s]'%base_dict[base]
	for base in PAM[::-1]:
		reverse_word+='[%s]'%Fasta_reverse(base_dict[base])
	if PAM_end==3:
		word=('.{%s}'%spacer_length)+word
		reverse_word=reverse_word+'.{%s}'%spacer_length
	elif PAM_end==5:
		word=word+('.{%s}'%spacer_length)
		reverse_word='.{%s}'%spacer_length+reverse_word
	rre=re.compile(word)
	rrv=re.compile(reverse_word)
	return rre,rrv

def load_gene(file):
	#info('Load %s'%file)
	gene_dic={}
	with open(file) as ff:
		for line in ff:
			line=line.strip()
			if line.startswith('>'):
				name=line.replace('>','')
				if not gene_dic.get(name):
					gene_dic[name]=[]
			else:
				line=line.upper()
				gene_dic[name].append(line)
	return gene_dic

def detect_repeat(spacer_list,target_sequence):
	if (target_sequence not in spacer_list) & (Fasta_reverse(target_sequence) not in spacer_list):
		spacer_list.append(target_sequence)
		spacer_list.append(Fasta_reverse(target_sequence))
		return spacer_list,True
	else:
		return spacer_list,False

def get_spacer(sss,PAM_end,spacer_length):
	if PAM_end==3:
		spacer=sss[:spacer_length]
	elif PAM_end==5:
		spacer=sss[len(PAM):]
	return spacer

def create_dict(mismatch):
	count_dict={}
	for i in range(mismatch+1):
		count_dict[str(i)]=0
	return count_dict

def get_info(mismatch_dic,spacer,x,y):
	mismatch_list=list(mismatch_dic[x].values())
	GC_content=count_GC(spacer)
	info_word='\t'.join(list(y)[:-1])
	loc=list(y)[-1]
	mismatch_word='\t'.join([str(i) for i in mismatch_list])
	total=sum(mismatch_list)-mismatch_list[0]
	return GC_content,info_word,mismatch_word,total,loc

def Scan(rre,rrv,gene_dic,PAM_end,spacer_length,base,window_size,window_start):
	spacer_list=[]
	count=0
	out=open('Temp.txt','w')
	fa_out=open('Temp.fa','w',encoding='utf-8')
	for key,value in gene_dic.items():
		value=''.join(value)
		sgRNA_f=rre.finditer(value,overlapped=True)
		sgRNA_s=rrv.finditer(value,overlapped=True)
		for i in sgRNA_f:			
			sss=i.group()
			spacer=get_spacer(sss,PAM_end,spacer_length)
			flag,loc=tell_window(spacer,base,window_size,window_start)
			if flag:
				if 'TTTT' not in spacer:
					spacer_list,flag=detect_repeat(spacer_list,spacer)
					if flag:
						count+=1
						start=i.start()
						end=i.end()
						name=key+'_S_'+str(count)
						out.write('%s\t%s\t%s\t%s\t%s\n'%(name,start,end,sss,loc))
						fa_out.write('>%s\n%s\n'%(name,spacer))
		count=0
		for i in sgRNA_s:			
			sss=Fasta_reverse(i.group())
			spacer=get_spacer(sss,PAM_end,spacer_length)
			flag,loc=tell_window(spacer,base,window_size,window_start)
			if flag:
				if 'TTTT' not in spacer:
					spacer_list,flag=detect_repeat(spacer_list,spacer)
					if flag:
						count+=1
						start=i.start()
						end=i.end()
						name=key+'_A_'+str(count)
						out.write('%s\t%s\t%s\t%s\t%s\n'%(name,start,end,sss,loc))	
						fa_out.write('>%s\n%s\n'%(name,spacer))
	out.close()
	fa_out.close()

def Cal_off(genome,result,mismatch,PAM_end,spacer_length):
	target='Temp.fa'
	if tell_plat():
		cmd='seqmap.exe %s %s %s Temp_map_out.txt /output_all_matches' % (mismatch,target,genome)
	else:
		cmd='./seqmap %s %s %s Temp_map_out.txt /output_all_matches' % (mismatch,target,genome)
	print(cmd)
	os.system(cmd)
	mismatch_dic={}
	with open('Temp.txt') as ff:
		for line in ff:
			sp=line.strip().split('\t')
			sg=sp[0]
			mismatch_dic[sg]=create_dict(mismatch)
	with open('Temp_map_out.txt') as ff:
		for line in ff:
			break
		for line in ff:
			sp=line.split('\t')
			sg=sp[3]
			m=sp[-2]
			if not mismatch_dic.get(sg):
				mismatch_dic[sg]=create_dict(mismatch)
			mismatch_dic[sg][m]+=1
	info_dic={}
	with open('Temp.txt') as ff:
		for line in ff:
			sp=line.strip().split('\t')
			info_dic[sp[0]]=sp
	count_list=list(list(mismatch_dic.values())[0].keys())
	with open(result,'w') as ot:
		if PAM_end==5:
			ot.write('SgID\tStart\tEnd\tSg_seq\tPAM\tSpacer\tLocation_of_target_base_in_spacer\tGC_content\tM'+'\tM'.join(count_list)+'\tTotal\n')
			for sgID,sp in info_dic.items():
				PAM_seq=sp[3][:len(PAM)]
				spacer=sp[3][len(PAM):]
				GC_content,info_word,mismatch_word,total,loc=get_info(mismatch_dic,spacer,sgID,sp)
				outword='%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (info_word,PAM_seq,spacer,loc,GC_content,mismatch_word,total)
				ot.write(outword)
		elif PAM_end==3:
			ot.write('SgID\tStart\tEnd\tSg_seq\tSpacer\tPAM\tLocation_of_target_base_in_spacer\tGC_content\tM'+'\tM'.join(count_list)+'\tTotal\n')
			for sgID,sp in info_dic.items():
				spacer=sp[3][:spacer_length]
				PAM_seq=sp[3][spacer_length:]				
				GC_content,info_word,mismatch_word,total,loc=get_info(mismatch_dic,spacer,sgID,sp)
				outword='%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (info_word,spacer,PAM_seq,loc,GC_content,mismatch_word,total)
				ot.write(outword)


def main():
	print(get_header())
	parser = argparse.ArgumentParser(description='CRISPRBEDesign Parameters')
	parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
	parser.add_argument('-i', '--input', type=str,  help='gene fast file', default='')
	parser.add_argument('-g', '--genome', type=str,  help='genome fast file', default='')
	parser.add_argument('-m', '--mismatch', type=int,  help='mismatch of sgRNA off-target prediction(default: 1)', default='1')
	parser.add_argument('-p', '--pam', type=str, help='Specify the pam to use for off-target analysis(default: NGG).', default='NGG')
	parser.add_argument('-s', '--spacer', type=int, help='Specify the spacer length.', default=20)
	parser.add_argument('-e', '--end', type=int, help='Specify the pam end(default: 3).', default=3)
	parser.add_argument('-W', '--Wsize', type=int, help='Set the base editing window size(default: 5).', default=5)
	parser.add_argument('-w', '--wstart', type=int, help='Set the start of base editing window(default: 4).', default=4)
	parser.add_argument('-t', '--target', type=str, help='Set the target base(default: A).', default='A')
	parser.add_argument('-o', '--output',  help='Output path', default='./result.txt')
	args = parser.parse_args()
	sys.stdout.flush()

	if not args.genome or not args.input:
		parser.print_help()
		exit(1)
		
	if os.path.exists(args.output):
		warn('%s exists,removing'%args.output)
		os.remove(args.output)	
			
	file=args.input
	genome=args.genome
	mismatch=args.mismatch
	PAM=args.pam
	PAM_end=args.end
	spacer_length=args.spacer
	result=args.output
	base=args.target
	window_start=args.wstart
	window_size=args.Wsize
	rre,rrv=mode(PAM_end,PAM,spacer_length)

	info('Loading gene file [%s]'%file)
	gene_dic=load_gene(file)
	info('Scaning gene file to get target site')
	Scan(rre,rrv,gene_dic,PAM_end,spacer_length,base,window_size,window_start)
	info('Caculating off-target effect.Please wait.')
	Cal_off(genome,result,mismatch,PAM_end,spacer_length)
	info('Caculating done! Removing temporary files.')
	os.remove('Temp.txt')
	os.remove('Temp.fa')
	os.remove('Temp_map_out.txt')
	info('Task done!')



if __name__ == '__main__':
	main()
