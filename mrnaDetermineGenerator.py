#!/usr/bin/python

from expressLevels import expressionLevel
from mrnalibrary import mrnaLibrary
from discrete_uniform_generator import discreteUniform
import numpy as np
import pdb
from mrnaGenerator import mrnaGenerator
# given a mrna template library and the target total library size
# this class can give the next transcript to be sequenced
class mrnaDetermineGenerator(mrnaGenerator):
	def __init__(self,mrnafile,poolSize,seed = 2147483647):
		mrnaGenerator.__init__(self,mrnafile,poolSize,seed)
		self.cur = 0
	def __iter__(self):
		return self
	# give the next ful length mrna for fragmentation
	def next(self):
		if self.cur < self.poolSz:
			tmp = self.cur
			self.cur += 1
			return self.mrnalib[self.expressLevel[tmp]]
		else:
			raise StopIteration()	
if __name__ == '__main__':
	mrnafile = 'mrna.txt'
	poolsize = 20
	mrnagen = mrnaDetermineGenerator(mrnafile,poolsize)
	for m in mrnagen:
		print(m)
        
    
