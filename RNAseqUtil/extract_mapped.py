#! /usr/bin/python

import sys
import os
from util.command.subproc import *
#from annotation import countMapped
def  source_parsing(file):
	dir = ''
	targets=[]
	genome = ''
	outputdir = ''
	for line in file:
		line = line.rstrip('\n')
		tuple = line.split(':')
		if tuple[0] =='directory':
			dir=tuple[1]
		elif tuple[0] == 'target':
			targets = targets + tuple[1:]
		elif tuple[0] == 'genome':
			genome = tuple[1]
		elif tuple[0] == 'outputdir':
			outputdir = tuple[1]
		else:
			print('undefined name' + tuple[0])
			sys.exit()
	return {'dir':dir,'targets':targets,'genome':genome,'outputdir':outputdir}
def extractRegionMapping(bamfile,region,outfile):
	command = 'samtools view ' +  bamfile +' '  + region + ' > ' + outfile 	
	os.popen(command)	
def extractSamRegion(samfile,region,outfile,rdlen = 50):
	handler = open(samfile,'r')
	regions={}
	for r in region:
		mrna,range = r.split(':')
		start,end =map(int,range.split('-'))
		len = end - start + 1
		start = start - rdlen + 1
		end = end 
		regions[mrna] = (start,end)
	outh = open(outfile,'w')
	for line in handler:
		if line[0] == '@':
			continue
		else:
			rec = line.split()
			if rec[2] in regions and int(rec[3]) >= regions[rec[2]][0] and int(rec[3]) <= regions[rec[2]][1]:
				#outh.write(rec[0]+'\t' + rec[2] + '\t' + rec[1] + '\t' + rec[3] + '\t' + rec[9] + '\t' + rec[5] + '\n')
				outh.write(line)
	outh.close()	
if __name__ == '__main__':
	#dir = '/home/stc/Documents/GenomicData/RNAseq/WBahou/THR102'
	from multiprocessing import Process	
	#bamfile = 'THR102_accepted_hits.bam'
	directory = '/home/luzhao/THR104_analysis/mapped'
	sample = ['A019','A023','N084','THR039','THR053']
	mutsample = ['THR100','THR101','THR102','THR104','THR116','THR136','THR137']
	#mrna = {'NM_203458ex5':['NM_203458ex5:1004-4962'],'NM_024408ex34':['NM_024408ex34:1057-11466','NM_024408ex34:49-51']}
	mrna = {'NM_001200001ex22':['NM_001200001ex22:1057-11466'],'NM_024408ex34':['NM_024408ex34:1057-11466']}
	outdir= '/home/luzhao/THR104_analysis/notch2/'
	procs = []
	for rna,r in mrna.items():
		print(rna,r)
		for s in sample + mutsample:
			p = Process(target = extractSamRegion,args=(os.path.join(directory,s + '.sam'),r,os.path.join(outdir,rna + '_' + s + '.sam')))
			p.start()
			procs.append(p)
		for p in procs:
			p.join()
	
