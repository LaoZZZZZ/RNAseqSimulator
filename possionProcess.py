#!/usr/bin/python

import os
import numpy as np
import scipy as cp
import math
# simulate a possition process
# the expected length of the fragment is the lambd
class PossionProcess(object):
	@staticmethod
	def cutPositions(lambd,end):
		pos = []
		start = 0  
		while True:
			tmp = PossionProcess.round(np.random.exponential(lambd))
			tmp_start = tmp + start
			tmp += start
			if tmp < end:
				pos.append(tmp)
				start = tmp_start	
			else:
				break
		pos.append(end)
		return pos
	@staticmethod
	def round(v):
		f = math.floor(v)
		if v - f > 0.5:
			return f + 1 
		else:
			return f 	

if __name__ == '__main__':
	np.random.seed(23)
	lambd = 200 
	cut = PossionProcess.cutPositions(lambd,2000)
	start = 1	
	frag = []
	while start < len(cut):
		frag.append(cut[start] - cut[start-1])
		start += 1
	print(cut,frag)		
 

