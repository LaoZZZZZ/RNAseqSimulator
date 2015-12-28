#!/usr/bin/python
import sys
import os
from alignment import alignment
from sampleInfo import sampleInfo
def sortSupp(file):
	try:
		handler = open(file,'r')
		supp = []
		for line in handler:
			rec = line.split()
			al = alignment(rec[0],rec[1],rec[2],rec[3],rec[4],rec[5])
			supp.append(al)
		handler.close()
		handler = open(file,'w')
		for e in sorted(supp):
			#print(e.__str__())
			handler.write(e.__str__())
		handler.close()	
	except Exception as err:
		print(err)
		sys.exit()
def sortAll(directory):
	dataset = sampleInfo()
	samples = dataset.samples()
	listSamples = []
	for k,v in samples.items():
		listSamples += v	
	for e in listSamples:
		os.chdir(os.path.join(directory,e))
		for supp in glob.glob('*.supp'):
			sortSupp(supp)	
if __name__ == '__main__':
	import glob
	directory = '/home/luzhao/THR104_analysis/mutationSam'
	sortAll(directory)
