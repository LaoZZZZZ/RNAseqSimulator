#!/usr/bin/python
from reads import read
from strMatch import findExactMismatches
class align_candidate(object):
	def __init__(self,rd):
		self.__rd = rd
		self.__match = []
	def getID(self):
		return self.__rd.getid()
	def getRead(self):
		return self.__rd
	def addmatch(self,m):
		self.__match = m
	def getMatch(self):
		return self.__match
	def __len__(self):
		return len(self.__match)
	def __str__(self):
		res = self.__rd.getid() + ' has %s matches : \n'%(len(self.__match))
		for head,tail in self.__match:
			res = res + ' head: %s \n tail: %s \n'%(head,tail)
		return res
	def __rec__(self):
		res = ''
		for head,tail in self.__match:
			res = res + str('%s%s'%(head.__str__(),tail.__str__())).rstrip('\n')
		return res
	def __eq__(self,other):
		ali1 = set()
		for head,tail in self.__match:
			ali1.add(head)
			ali1.add(tail)
		ali2 = set()
		for head,tail in other.getMatch():
			ali2.add(head)
			ali2.add(tail)
		return not (ali1 ^ ali2)
	def __hash__(self):
		h = 0
		for head,tail in self.__match:
			h = h + head.__hash__() + tail.__hash__()
		return h
		#return findExactMismatches(self.__rd.getseq(),other.getRead().getseq()) <= 3	

