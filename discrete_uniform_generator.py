#!/usr/bin/python

import numpy as np

# give a list of random variates between the [low,high] interval
# 
class discreteUniform(object):
	@staticmethod
	def setSeed(seed):
		np.random.seed(seed)
	@staticmethod
	def next(low,high = None,size=None):
		return np.random.random_integers(low,high,size)


if __name__ == '__main__':
	for i in range(10):
		print(discreteUniform.next(2))
