#!/usr/bin/python
import os
import sys

from reads import read 

from read_func import *
from scan_reads import loadSam
from scan_reads import singleSample as scan
import gzip
dir = '/home/luzhao/calr_analysis'
readf = open('/home/luzhao/calr_analysis/unmapped/THR104.fq','r')
dict = loadSam(dir,'samfile/THR104')
#scan('/home/luzhao/calr_analysis/script',dir,'short','/home/luzhao/calr_analysis/script',20)	
total = 0
totalalign = 0 
record = []
totalreads = 0.0
totalunmapped = 0.0
for line in readf:
	total = total + 1
	record.append(line.rstrip('\n'))
	if total % 4 == 0:
		totalreads = totalreads + 1	
		rd = read(record[0],record[1],record[3])
		if not rd.getid() in dict:	
			totalunmapped = totalunmapped + 1
		else:
			totalalign = totalalign + 1
		record = []

print('Total reads: %s\nUnmapped:%s\nmapped:%s\n'%(totalreads,totalunmapped,totalalign))		



