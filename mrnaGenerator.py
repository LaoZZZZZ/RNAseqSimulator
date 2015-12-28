#!/usr/bin/python

from expressLevels import expressionLevel
from mrnalibrary import mrnaLibrary
from discrete_uniform_generator import discreteUniform
import numpy as np
import pdb
import json
# given a mrna template library and the target total library size
# this class can give the next transcript to be sequenced
class mrnaGenerator(object):
	def __init__(self,mrnafile,poolSize,seed = 2147483647):
		#pdb.set_trace()
		#print(mrnafile)
		self.mrnalib = mrnaLibrary(mrnafile)
		self.poolSz = poolSize
		self.seed = seed
		np.random.seed(self.seed)
		print('In the process of library preparation')
		self.expressLevel = expressionLevel(len(self.mrnalib),self.poolSz)
		print('The size of the transcriptome is %s'%(self.poolSz))
		print('Library preparation finished')
	def __iter__(self):
		return self
	# give the next ful length mrna for fragmentation
	def next(self):
		raise StopIteration()	
	# get the copy numbers for each mrna transcripts in dictionary data structure				
	def getExpressionlevels(self):
		levels = self.expressLevel.expressLevels()
		res = {}
		absolute = 0.0
		totalbp = 0.0
		for i in levels.keys():
			length = len(self.mrnalib[i][1])
			copy = levels[i]
			relative = length*copy
			totalbp += relative	
			res[self.mrnalib[i][0]] = [length,copy,relative]
		for k in res.keys():
			res[k][1] /= float(self.poolSz)
			res[k][2] /= totalbp  
		return res
	def __len__(self):
		return self.poolSz 
	def dump(self,outfile):
		expressions = self.getExpressionlevels()
		absolute = 0.0
		relative = {}
		for k,v in expressions.items():
			relative[k] = v[0] * v[1]
			absolute += relative[k]
		handler = open(outfile,'w')
		handler.write('mrna,length,absolute expression,relative expression\n')
		for k in expressions.keys():
			handler.write('%s,%s,%s,%s\n'%(k,str(expressions[k][0]),str(expressions[k][1]/float(self.poolSz)),str(relative[k]/absolute)))
		handler.close()
if __name__ == '__main__':
	mrnafile = 'mrna.txt'
	poolsize = 20
	mrnagen = mrnaGenerator(mrnafile,poolsize)
	for m in mrnagen:
		print(m)
        
    
