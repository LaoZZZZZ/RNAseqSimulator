#!/usr/bin/python

import os
import sys
import gzip
import math
from extract_mapped import source_parsing
from sampleInfo import sampleInfo
def ErrorProb(c):
	return math.pow(10,-1*ord(c)/10)
def qualityLog(log,qualities):
	for i in range(len(qualities)):
		log[i][qualities[i]] += 1
	return log
def printLog(log,file):
	localfile = open(file,'w')
	index = 1	
	header = 'offset'
	for j in range(94):
		header += ',%s'%(j)
	header += '\n'
	localfile.write(header)	
	for e in log:
		line = '%s,'%(index)
		index += 1
		if sum(e.values()) == 0:
			continue
		for k,v in e.items():
			line += str(v) 
			line += ','
		line += '\n'
		localfile.write(line)
	localfile.close()
def qualityStatistics(sample,outdir):
	samples = sampleInfo()
	log = []
	for i in range(55):
		log.append(dict())
		for j in range(94):
			log[i][j+33] = 0
	for e in samples.readsData(sample):	
		dest = os.path.join(task['dir'],sample + '/' +e)
		fhand = gzip.open(dest,'r')
		total = 1 
		for line in fhand:
			if not total %  4:
				qualities = map(ord,line.rstrip())
				log = qualityLog(log,qualities)		
			total = total + 1
	output = sample + '.csv'
	printLog(log,outdir + '/' + output) 
if __name__ == '__main__':
	task = source_parsing(sys.stdin)
	for e in task['targets']:
		qualityStatistics(e,task['outputdir'])
