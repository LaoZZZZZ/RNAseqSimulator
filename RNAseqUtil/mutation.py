#!/usr/bin/python
from align_candidate import align_candidate
from refGeneStruct import refGene
import hashlib
import binascii
class mutation(object):
	def __init__(self,mutation_type,chrom,pos,mutSeq,matches):
		# 1 is insertion, 2 is deletion 
		self.__mutation_type = int(mutation_type)
		self.__mutChr = chrom
		self.__mutPos = int(pos)
		self.__mutSeq = mutSeq
		self.__matches = matches
		self.__normalmatches = []
		self.ids = set()
		for i in self.__matches:
			self.ids.add(i.getID())
	def mutationType(self):
		return self.__mutation_type
	def mutationPos(self):
		return self.__mutPos
	def mutationStartPos(self,refgene):
		if refgene.orient() == '+':
			return self.__mutPos
		else:
			if self.__mutation_type == 1 or self.__mutation_type == 3:
				return self.__mutPos
			else:
				#print(self.__mutPos + len(self.__mutSeq) - 1)
				return self.__mutPos + len(self.__mutSeq) - 1
	def mutationEndPos(self,refgene):
		if refgene.orient() == '+':
			if self.__mutation_type == 1 or self.__mutation_type == 3:
				return self.__mutPos  
			else:
				return self.__mutPos + len(self.__mutSeq)  - 1
		else:
			return self.__mutPos
	def setChrStart(self,start):
		self.__chrStart = start
	def setChrEnd(self,end):
		self.__chrEnd = end
	# return a close interval[start,end]
	def setChrMutRange(self,refgene):
		self.__chrStart = refgene.getChrPos(self.mutationStartPos(refgene))
		self.__chrEnd = refgene.getChrPos(self.mutationEndPos(refgene))
	def getChrMutRange(self):
		return (self.__chrStart,self.__chrEnd)
	def mutationChr(self):
		return self.__mutChr
	def setRna(self,mrna):
		self.__mutChr = mrna
	def setMutationPos(self,pos):
		self.__mutPos = pos
	def mutationSeq(self):
		return self.__mutSeq
	def supportReads(self):
		rd = []
		for e in self.__matches:
			rd.append(e.getRead())
		return rd
	def NumOfSupportReads(self):
		return len(set((self.__matches)))
	def getSupport(self):
		return self.__matches
	def getNormalSupport(self):
		return self.__normalmatches
	def addNormalSupport(self,match):
		if type(match) is list:
			for ind in match:
				self.__normalmatches.append(ind) 
	        elif isinstance(match,align_candidate):
			if not match.getID() in ids:
				self.__normalmatches.append(match)
		else:
			print(type(match))
			raise Exception('invalid type\n')
	def filterDuplicateReads(self,matches):
		ids = set(map(align_candidate.getID,matches))
		tmp = []
		for i in matches:
			if i.getID() in ids:
				tmp.append(i)
				ids.remove(i.getID())
		return tmp
	def addSupport(self,match):
		if type(match) is list:
			for ind in match:
				self.__matches.append(ind) 
	        elif isinstance(match,align_candidate):
			if not match.getID() in ids:
				self.__matches.append(match)
		else:
			print(type(match))
			raise Exception('invalid type\n')
	def setChromosome(self,chr):
		self.__chromo = chr
	def getChromo(self):
		return self.__chromo
	def setChrPos(self,pos):
		self.__chrPos = pos
	def getChrPos(self):
		return self.__chrPos
	def setGene(self,gene):
		self.__gene = gene
	def getGene(self):
		return self.__gene
	#def setUniqMapp(self,num):
	#	self.__uniqueCov = num
	#def setOverallMapp(self,total):
	#	self.__totalCov = total
	def __len__(self):
		return len(self.__mutSeq)
	def __str__(self):
		self.__matches = self.filterDuplicateReads(self.__matches)
		self.__normalmatches = self.filterDuplicateReads(self.__normalmatches)
		res = ''
		res = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,'%(self.__mutation_type,self.__mutChr,self.__mutPos,
						self.__chromo,self.__chrStart,self.__chrEnd,self.__gene,len(self.__mutSeq),
						self.__mutSeq,len(set(self.__matches)),len(self.__matches),
						len(set(self.__normalmatches)),len(self.__normalmatches),
						float(len(set(self.__matches)))/(float(len(set(self.__normalmatches)))+ float(len(set(self.__matches)))),
						float(len(self.__matches))/(float(len(self.__matches)) + float(len(self.__normalmatches))))
		for e in self.__matches:
			res = res + e.getID() + ' ' 
		res = res + ','
		for e in self.__normalmatches:
			res = res + e.getID() + ' '
		return res 
	def __hash__(self):
		return int(hashlib.md5(binascii.a2b_qp('%s%s%s'%(self.__mutChr,self.__mutPos,self.__mutSeq))).hexdigest(),16)
	def __lt__(self,other):
		try:
			if type(self) != type(other):
				raise Exception('can not compare with other type(%s) of instance with mutation class instance'%(type(other)))
			return self.NumOfSupportReads() < other.NumOfSupportReads()
		except Exeption as err:
			print(err)
			import sys
			sys.exit(0)
	
	def __gt__(self,other):
		try:
			if type(self) != type(other):
				raise Exception('can not compare with other type(%s) of instance with mutation class instance'%(type(other)))
			return self.NumOfSupportReads() > other.NumOfSupportReads()
		except Exeption as err:
			print(err)
			import sys
			sys.exit(0)
	def __eq__(self,other):
		if type(self) != type(other):
			return False	
		return self.__hash__() == other.__hash__()
