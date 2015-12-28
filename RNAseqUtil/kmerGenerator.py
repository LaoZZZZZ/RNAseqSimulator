#!/usr/bin/python

import os
import sys
from reads import read
class kmerGenerator:
	def __init__(self,klen,maxN = 2,Nstrat = 'A'):
		if klen %2 == 0:
			print('only the length of kmer be able to odd!\n')
			sys.exit(0)
		self.__klen = klen
		self.__dict = {'A':0,'a':0,'C':3,'c':3,'G':2,'g':2,'T':1,'t':1}
		self.__maxN = maxN
		#self.__dict = {'A':0,'a':0,'C':1,'c':1,'G':2,'g':2,'T':3,'t':3}
		self.__Nstrat = Nstrat
		if Nstrat.upper() in ['A','T']:
			self.__Ncomp = 1 - self.__dict[Nstrat]
		else:
			self.__Ncomp = 5 - self.__dict[Nstrat]
	def bitwise(self,seq,start):
		fw = 0
		rv = 0
		itr = 0
		while itr < self.__klen:
			c = seq[start+itr]
			fw = fw << 2
			if c in self.__dict:
				fw |= self.__dict[c]
			elif c in ['n','N']:
				fw |= self.__dict[self.__Nstrat]
			else:
				raise Exception('There are some other character other than A,C,G,T,N which is ' +c)
			rv = rv <<2
			c = seq[start+self.__klen-itr-1]
			if c in self.__dict:
				if c.upper() in ['A','T']:
					rv |= ( 1 - self.__dict[c])
				else:
					rv |= (5 - self.__dict[c])
				#rv |= (3 - self.__dict[c])
			elif c in ['n','N']:
				#rv |= (3 - self.__dict[self.__Nstrat])
				rv |= self.__Ncomp
			else:
				raise Exception('There are some other character other than A,C,G,T,N which is ' +c)
			itr = itr + 1
		return [fw,rv]
			
	def generateKmer(self,rd):
		tmp = {}
		if len(rd) >= self.__klen and (rd.getseq().count('n') + rd.getseq().count('N') <= self.__maxN):
			seq = rd.getseq()
			for itr in range(len(seq) - self.__klen + 1):
				res = self.bitwise(seq,itr)
				if res[0] in tmp:
					tmp[res[0]] = tmp[res[0]] +  1
				else:
					tmp[res[0]] = 1
				if res[1] in tmp:
					tmp[res[1]] = tmp[res[1]] +  1
				else:
					tmp[res[1]] = 1
		return tmp	
