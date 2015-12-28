#!/usr/bin/python

from expressLevels import expressionLevel
from mrnalibrary import mrnaLibrary
from discrete_uniform_generator import discreteUniform
import numpy as np
import pdb
from mrnaGenerator import mrnaGenerator
# given a mrna template library and the target total library size
# this class can give the next transcript to be sequenced
class mrnaRandomGenerator(mrnaGenerator):
	def __init__(self,mrnafile,poolSize,seed = 2147483647):
		mrnaGenerator.__init__(self,mrnafile,poolSize,seed)
	def __iter__(self):
		return self
	# randomly give the next ful length mrna for fragmentation
	def next(self):
		ind = discreteUniform.next(0,self.poolSz - 1,1)
		return self.mrnalib[self.expressLevel[ind]]
if __name__ == '__main__':
	mrnafile = 'mrna.txt'
	poolsize = 20
	mrnagen = mrnaRandomGenerator(mrnafile,poolsize)
	print(len(mrnagen))
	#for m in mrnagen:
	#	print(m)
        
    
