#!/usr/bin/python

import os
import sys
from fragment import Fragment
# filter those fragment that does not satisfy the condition
# currently, we only keep those fragment whose length is between [lower,higher]
class filterFragment(object):
	def __init__(self,lower,higher):
		self.lower = lower
		self.higher = higher
	# input: a list of fragment, 
	def filterFragment(self,f):
		return len(f) >= self.lower and len(f) <= self.higher
	def filterFragments(self,fs):
		return [ e for e in fs if self.filterFragment(e)] 


if __name__ == '__main__':
	f1 = Fragment('chr1',0,20)		
	f2 = Fragment('chr2',20,100)
	fs = [f1,f2,Fragment('chr3',100,10)]
	filterF = filterFragment(40,100)
	print(filterF.filterFragment(f1),filterF.filterFragment(f2))
	print(filterF.filterFragments(fs))
