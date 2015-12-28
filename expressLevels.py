#!/usr/bin/python
from mrnalibrary import mrnaLibrary as mrnaLib
from discrete_uniform_generator import discreteUniform as disUnif


import os
import numpy as np

# given a library size and the total number of random sample sz
# give the expression level for each mrna labeled [0 - libsz-1]
class expressionLevel(object):
	def __init__(self,libsz,total):
		self.libsz = libsz
		self.total = total
		self.generateSettings()
	def generateSettings(self):
		self.pool = disUnif.next(0,self.libsz-1,self.total)
		self.express = {}
		for ind in self.pool:
			self.express.setdefault(ind,0)
			self.express[ind] += 1
	def expressLevels(self):
		return self.express
	def getPool(self):
		return self.pool
	def __getitem__(self,ind):
		return self.pool[ind]

if __name__ == '__main__':
	sz = 10
	total = 500
	level = expressionLevel(sz,total)
	express = level.expressLevels()
	print(express)
	print(level[10:14])
    
  
