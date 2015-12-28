#!/usr/bin/python

class Sam(object):
	def __init__(self,rdname,pos,refname,cigar = '',mapq = 0,flag = 0):
		self.qname = rdname
		self.pos = pos
		self.refname = refname
		self.cigar = cigar
		self.mapq = mapq
		self.flag = flag
		self.nref = '*'
		self.npos = 0
		self.tlen = 0
		self.seq = ''
		self.qua = ''
	def __repr__(self):
		return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(self.qname,self.flag,self.refname,self.pos,self.mapq,self.cigar,self.nref,self.npos,self.tlen,self.seq,self.qua)
	def __str__(self):
		return self.__repr__() 
	def setCigar(self,cigar):
		self.cigar = cigar
	def setSeq(self,seq):
		self.seq = seq
	def setQua(self,qua):
		self.qua = qua
	def setMapQ(self,mapQ):
		self.mapq = mapQ
	def setFlag(self,f):
		self.flag = f


	
