#!/usr/bin/python
import numpy as np
import os
from scipy.stats import bernoulli
from scipy.stats import norm
from discrete_uniform_generator import discreteUniform
import math
#given the prob of sequence error
# generate the base pair based on the offset and the original base pair info
class sequenceError(object):
	# currently we assume the quality follows a norm distribution with mean = 90,std= 10
	# the position effect has not be considered which should be incoporated later for better simulation
	def __init__(self,locale = 85,scale = 15):
		self.norm = norm(locale,scale)
		self.set = set(['A','C','G','T'])
	#the error probability is given by 10^(-1*qua/10)
	# then perform a bernoulli experiment with the success prob = error probability to decide wheather to give sequence error	
	def basePair(self,bp,offset):
		qua = self.norm.rvs(1)
		if qua > 123:
			qua = 123
		elif qua < 66:
			qua =66
		t = qua-34
		if t < 0:
			t = 0
		r = bernoulli.rvs(math.pow(10,-1*t/10),size=1)
		if r:
			choice = np.random.choice(list(self.set-set(bp)))	
			return (choice,int(qua))
		else:
			return (bp,int(qua))


if __name__ == '__main__':

	seq = sequenceError(90,5)
	total = 1000000
	err = 0.0
	for i in range(total):
		bp = seq.basePair('A',10)[0]
		if bp != 'A':
			err += 1 
 	print(err,err/total )
