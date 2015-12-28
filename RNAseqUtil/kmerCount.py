#!/usr/bin/python

import os
import sys
import struct

from reads import read
from kmerGenerator import kmerGenerator as generator

def update(dict,extra):
	for k,v in extra.items():
		if k in dict:
			if dict[k] == 255 or dict[k] + v > 255:
				dict[k] = 255
			else:
				dict[k] = dict[k] + v
		else:
			if v > 255:
				dict[k] = 255
			else:
				dict[k] = v
	return dict
def write(dict,output):
	fhandler = open(output,'wb')
	for k,v in dict.items():
		fhandler.write(struct.pack('<qB',k,v))
	fhandler.close()
def kmerCount(inputdir,target,klen,outputdir):
	gen = generator(klen)
	input = os.path.join(inputdir,target+'.fq')
	input = open(input,'r')
	output = os.path.join(outputdir,target + '.kmer')
	total = 0
	lines = []
	kmers = {}
	for line in input:
		total = total + 1
		lines.append(line.rstrip('\n'))
		if total % 4 == 0:
			if lines[2] != '+':
				print('invalid fq file')
				sys.exit(0)
			rd = read(lines[0],lines[1],lines[3])
			kmers = update(kmers,gen.generateKmer(rd))
			lines = []
	write(kmers,output) 	

if __name__ == '__main__':
	inputdir = '/home/luzhao/calr_analysis/kmer'
	target = ['A019','A023','THR039','THR053','THR100','THR101','THR102','THR103','THR116','THR136','THR137']
	klen = 23
	outputdir = '/home/luzhao/calr_analysis/kmer'
	for e in target:
		kmerCount(inputdir,e,klen,outputdir)
