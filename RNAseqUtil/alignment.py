#!/usr/bin/python
import hashlib
import binascii
import re
class alignment(object):
	def __init__(self,id,chrom,orient,pos,seq,score = 0,cigar = ''):
		self.__id = id
		self.__chr = chrom
		self.__dir = orient
		self.__pos = pos
		self.__seq = seq.upper()
		self.__cigar = cigar
		self.__score = score
	def getID(self):
		return self.__id
	def getSeq(self):
		return self.__seq
	def getOrient(self):
		return int(self.__dir)
	def getPos(self):
		return int(self.__pos)
	def getChr(self):
		return self.__chr
	def getCigar(self):
		return self.__cigar
	def setPos(self,pos):
		self.__pos = pos 
	def setOrient(self,orient):
		self.__dir = orient
	def setChr(self,chr):
		self.__chr = chr
	def getScore(self):
		return self.__score
	def setCigar(self,newcigar):
		self.__cigar = newcigar	
	def __str__(self):
		return '%s\t%s\t%s\t%s\t%s\t%s\n'%(self.__id,self.__chr,self.__dir,self.__pos,self.__seq,self.__cigar)
	def __eq__(self,other):
		if type(other) != type(self):
			return False
		else:
			#return self.__pos == other.getPos() and self.__dir == other.getOrient() and self.__chr == other.getChr()
			return self.__seq == other.getSeq()

	def __len__(self):
		return len(self.__seq)
	def __lt__(self,other):
		try:
			if self.__chr != other.getChr():
				raise Exception( 'can not compare two alignments from different reads')
			return  self.__pos < other.getPos()
		except Exception as err:
			print(err)
			import sys
			sys.exit(0)
        
	def __le__(self,other):
		try:
			if self.__id != other.getID():
				raise Exception( 'can not compare two alignments from different reads')
			return  self.__pos <= other.getPos()
		except Exception as err:
			print(err)
			import sys
			sys.exit(0)
	def __hash__(self):
		#return int(hashlib.md5(binascii.a2b_qp('%s%s%s'%(self.__chr,self.__dir,self.__pos))).hexdigest(),16) 
		return int(hashlib.md5(binascii.a2b_qp(self.__seq)).hexdigest(),16) 
	def getRange(self):
		start = self.__pos
		end = self.__pos
		cut = re.split('M|I|D',self.__cigar)[0:-1]
		off = 0
		roff = 0
		moff = 0
		for p in cut:
			off = off + len(p)
			if self.__cigar[off] == 'M':
				end += int(p)
				roff = roff + int(p)
				moff = moff + int(p)
				off = off + 1
			elif self.__cigar[off] == 'I':
				roff = roff + int(p)
				off = off + 1
			elif self.__cigar[off] == 'D':
				end += int(p)
				off = off + 1
				moff = moff + int(p)
		return (start,end-1)
	def getMutations(self):
		start = self.__pos
		end = self.__pos
		cut = re.split('M|I|D',self.__cigar)[0:-1]
		off = 0
		roff = 0
		moff = 0
		mutations = set() 
		for p in cut:
			off = off + len(p)
			if self.__cigar[off] == 'M':
				end += int(p)
				roff = roff + int(p)
				moff = moff + int(p)
				off = off + 1
			elif self.__cigar[off] == 'I':
				mutations.add((end,end,'I'))
				roff = roff + int(p)
				off = off + 1
			elif self.__cigar[off] == 'D':
				mutations.add((end,end+int(p)-1,'D'))	
				end += int(p)
				off = off + 1
				moff = moff + int(p)
		return mutations

