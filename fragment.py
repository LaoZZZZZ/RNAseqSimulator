#!/usr/bin/python

class Fragment(object):
	def __init__(self,transcriptID,start,length):
		self.id = transcriptID
		self.start = int(start)
		self.length = int(length)
	def __len__(self):
		return self.length
	def startPos(self):
		return self.start
	def __str__(self):
		return '%s:%s,%s'%(self.id,self.start,self.length)
	def __repr__(self):
		return '%s:%s,%s'%(self.id,self.start,self.length)
	def transcriptID(self):
		return self.id
