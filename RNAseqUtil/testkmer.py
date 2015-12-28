#!/usr/bin/python

import os
import sys
from reads import read
from kmerGenerator import kmerGenerator

if __name__ == '__main__':
	
	gen = kmerGenerator(23)
	total = 0
	rec = []
	for line in sys.stdin:
		total = total + 1
		rec.append(line.rstrip('\n'))
		if total % 4  == 0:
			rd = read(rec[0],rec[1],rec[3])
			if rec[2] 1= '+':
				print('invalid parsing of fq file')
				sys.exit(0)
					
	
