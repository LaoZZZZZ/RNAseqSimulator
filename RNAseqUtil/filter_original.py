#! /usr/bin/python
import gzip
import os
import sys
from multiprocessing import Process,Lock

from extract_mapped import source_parsing
from reads import read
from read_func import split, getFrontEnd
from scan_reads import loadSam
from scan_reads import trim
def filterOrigin(fqdir,samdir,target,output):
	inputfq = os.path.join(fqdir,target+'.fq')
	dict = loadSam(samdir,target)		
	outputfile = os.path.join(output,target+'.fq')
	outhandler = open(outputfile,'w')
	inhandler = open(inputfq,'r')
	total = 0
	rd = []
	total = 0
	for line in inhandler:
		total = total + 1
		rd.append(line.rstrip('\n'))
		if total %  4 == 0:
			if not '+' in rd[2]:
				print('incorrect fq  format')
				sys.exit(0)
			if trim(rd[0]) in dict: 
				rdbody = read(rd[0],rd[1],rd[3])
				outhandler.write(rdbody.__str__()+'\n')
			rd = []
	outhandler.close()
if __name__ == '__main__':
	task = source_parsing(sys.stdin)
	print(task)
	samdir = '/home/luzhao/calr_analysis/region_mapped' 
	output = '/home/luzhao/calr_analysis/kmer'
	if not os.path.isdir(output):
		os.system('mkdir ' + output)		
	for target in task['targets']:
		print(target)
		Process(target=filterOrigin,args=(task['dir'],samdir,target,output)).start()
	print('finished filtering')
