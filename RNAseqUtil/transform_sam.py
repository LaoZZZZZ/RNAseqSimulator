#! /usr/bin/python

import sys
import os
from extract_mapped import source_parsing
from util.command.subproc import *
if __name__ == '__main__':
	print('executing')
	data = source_parsing(sys.stdin)
	dir = data['dir']
	targets=data['targets']
	suffix = '.sam'
	output = '/home/luzhao/calr_analysis/samfile'
	if not os.path.isdir(output):
		os.system('mkdir ' + output)
	for e in targets:
		tmp = os.path.join(dir,e)
		for f in os.listdir(tmp):
			if os.path.isfile(os.path.join(tmp,f)) and f.endswith('bam'):
				input = os.path.join(tmp,f)
				outputfile = os.path.join(output,e+suffix)
				command = 'samtools view ' +  input + ' > ' + outputfile 	
				os.popen(command)	
