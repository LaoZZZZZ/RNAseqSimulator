#!/usr/bin/python
from alignment import alignment
import os
class suppFileIterator:
	def __init__(self,filename):
		self.__handler = open(filename,'r')
	def __iter__(self):
		return self
	def next(self):
		line = self.__handler.readline().rstrip()
		while line and '@' == line[0]:
			line = self.__handler.readline().rstrip()
		if line:
			rec = line.split()
			if len(rec) > 4 :
				al = alignment(rec[0],rec[1],rec[2],rec[3],rec[4],0,rec[5])
				return al
		else:
			raise StopIteration		
		
if __name__ == '__main__':
	file = '/home/luzhao/multipleHits/THR104_first_round_bestAlignments_test.sam'
	sam = suppFileIterator(file)
	for line in sam:
		print(line.__str__())
		break
	
