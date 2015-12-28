#!/usr/bin/python

from kmerCount import kmerCount
from kmerCount import write
from reads import read
from kmerGenerator import kmerGenerator as kgen
import os
import sys


def parse_config(file):
	chr = ''
	klen = 0
	out = ''
	try:
		for line in file:
			rec = line.rstrip('\n').split(':')
			if rec[0].lower() in ['genome','chr']:
				chr = rec[1]
			elif rec[0].lower() in ['len','klen','kmer']:
				klen = int(rec[1])
			elif rec[0].lower() in ['output','out']:
				out = rec[1]
			
	except Exception as error:
		print(error)
		sys.exit(0)
	return [chr,klen,out]
def readfa(file):
	try:
		chr = [] 
		body = [] 
		f = open(file,'r')
		tmp = ''
		for line in f:
			line = line.rstrip('\n')
			if line[0] == '>' and not tmp:
				tmp = line
				chr.append(tmp)
			elif line[0] == '>' and tmp:
				body.append('')
				tmp = ''
			elif line[0] != '>' and tmp :
				body.append(line)
				tmp = ''
			else:
				raise Exception(file)
		return {'chr':chr,'seq':body}
	except Exception as err:
		print(err)
				 	
def refKmer(input):
	try:
		chr = input[0]
		klen = int(input[1])
		out = input[2]
		gen = kgen(klen)
		seq = readfa(chr)
		dict = gen.generateKmer(read('1',seq['seq'][0],''))
		write(dict,out)
	except Exception as error:
		print(error)
		sys.exit(0)
if __name__ == '__main__':
	#input = parse_config(sys.stdin)
	#print(input)
	#refKmer(input)
	import os
	direc = '/home/luzhao/THR104_analysis/chromo/ref'
	for f in os.listdir(direc):
		readfa(os.path.join(direc,f))	
