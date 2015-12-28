#!/usr/bin/python
from alignment import alignment
import os
class SESamfileIterator:
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
			if len(rec) > 9:
				al = alignment(rec[0],rec[2],rec[1],rec[3],rec[9],rec[4],rec[5])
				return al
		else:
			raise StopIteration		
		
if __name__ == '__main__':
	import glob
	filename = '/home/luzhao/THR104_analysis/mapped/'
	for e in glob.glob(filename + '/' + 'NM*'):
		print(os.path.join(filename,e))
		sam = SESamfileIterator(os.path.join(filename,e))
		for line in sam:
			print(line.__str__())
			break
	
