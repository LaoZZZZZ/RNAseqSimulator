#!/usr/bin/python

import os
import sys
from sets import Set
from extract_mapped import source_parsing

def sort_bam(inputdir,target):
	sort = 'samtools sort '
	sam2bam='samtools view -bS '
	bam2sam = 'samtools view -h '
	samfile = os.path.join(inputdir,target+'.sam')
	bamfile = os.path.join(inputdir,target+'_sorted')
	output = os.path.join(inputdir,target+'_sorted.sam')
	tmp = os.path.join(inputdir,'tmp.bam')
	if os.path.isfile(samfile):
		command = sam2bam + samfile + ' > ' + tmp  
		os.popen(command)
		command = sort + ' ' + tmp  + ' ' + bamfile  
		os.popen(command)
		command = bam2sam + ' ' + bamfile  + '.bam ' +  '> ' + output
		os.popen(command)
		os.system('rm ' + tmp)
		  
if __name__ == '__main__':
	task = source_parsing(sys.stdin)
	print(task)
	inputdir = task['dir']
	for target in task['targets']:
		 sort_bam(inputdir,target) 
