#! /usr/bin/python
import gzip
import os
import sys
from multiprocessing import Process,Lock

from extract_mapped import source_parsing
from reads import read
from read_func import split, getFrontEnd
from scan_reads import loadSam
from filter_original import filterOrigin
if __name__ == '__main__':
	task = source_parsing(sys.stdin)
	print(task)
	samdir = '/home/luzhao/calr_analysis/mapped' 
	output = '/home/luzhao/calr_analysis/final_sam'
	if not os.path.isdir(output):
		os.system('mkdir ' + output)		
	for target in task['targets']:
		print(target)
		Process(target=filterOrigin,args=(task['dir'],samdir,target,output)).start()
	print('finished filtering')
