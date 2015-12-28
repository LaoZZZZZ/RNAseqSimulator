import math
import sys
from sets import Set
import numpy as np
import pdb
# class holds single reads    
class read(object):
	def __init__(self,id,seq,qua):
		self.__id = id
		self.__seq = seq
		self.__qua = qua
		if len(self.__seq) != len(self.__qua):
			print "the length of quality is not equal to the length of read sequence!"
			sys.exit(0)
	def __len__(self):
		return len(self.__seq)
	def seq(self):
		return self.__seq
	def name(self):
		return self.__id
	def qual(self):
		return self.__qua
	def empty(self):
		return not self.__seq
	def __repr__(self):
		return '@%s\n%s\n+\n%s'%(self.__id,self.__seq,self.__qua)
	def __str__(self):
		return self.__repr__()
	def complem(self):
		rc = ''
		for c in self.__seq.upper():
			if c == 'A':
				rc += 'T'
			elif c  == 'C':
				rc += 'G'
			elif c == 'T':
				rc += 'A'
			elif c == 'G':
				rc += 'C'
			elif c == 'N':
				rc += 'N'
			else:
				print('invalid squence')
				sys.exit()
		self.__seq = rc
	def reverseComp(self):
		rc = ''
		for c in reversed(self.__seq.upper()):
			if c == 'A':
				rc += 'T'
			elif c  == 'C':
				rc += 'G'
			elif c == 'T':
				rc += 'A'
			elif c == 'G':
				rc += 'C'
			elif c == 'N':
				rc += 'N'
			else:
				print('invalid squence')
				sys.exit()
		self.__seq = rc
				
