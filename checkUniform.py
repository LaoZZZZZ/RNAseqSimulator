#!/usr/bin/python

import sys
from drawQQplot import qqplotpdf
from rnaSimulator import drive
def getStartPos(file):
	res = []
	for line in open(file,'r'):
		rec = line.split()
		res.append(int(rec[3]))
	return res
def checkUniform(mrnafile,samfile,outfile,title):
	positions = getStartPos(samfile)
	mrna = open(mrnafile,'r').read().split('\n')[1].strip()
	length = len(mrna)
	print(length,len(positions))
	positions = [float(e)/length for e in positions]
	qqplotpdf(positions,'uniform',outfile,title)		
if  __name__ == '__main__':
	
	samfile = 'RNASimulation.sam'
	fold = int(sys.argv[1])
	drive(fold)
	title = '%s fold'%(str(fold))
	outfile = 'startPosQQ%s.pdf'%(str(fold))
	checkUniform('mrna.txt',samfile,outfile,title)
