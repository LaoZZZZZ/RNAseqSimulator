#!/usr/bin/python

# all position start from 0
class refGene(object):
	dictmap = {'A':'T','a':'T','T':'A','t':'A','G':'C','g':'C','C':'G','c':'G','N':'N','n':'N'}
	def __init__(self,mrna_id,chr,gstart,gend,cdstart,cdend,orient,exNum,exStarts,exEnds,altername):
		self.__mrnaid = str(mrna_id)
		self.__chr = chr
		self.__gstart = int(gstart )
		self.__gend = int(gend)
		self.__cdstart = int(cdstart)
		self.__cdend = int(cdend)
		self.__orient = orient
		self.__exNum = int(exNum)
		self.__exStarts = map(int,exStarts)
		self.__exEnds = map(int,exEnds)
		self.__alternativeName = altername
		# change the mrna id by addding the $ + number of exons
		self.checkSanity()
		self.__mrnaid = self.__mrnaid + 'ex' + str(self.__exNum) 
		self.calculatePos()
		self.codingRegion()
	def __len__(self):
		return self.__len
	def calculatePos(self):
		self.__exlen = []
		self.__len = 0
		for i in range(len(self.__exStarts)):
			self.__len += (self.__exEnds[i] - self.__exStarts[i])
			self.__exlen.append(self.__exEnds[i] - self.__exStarts[i])	
	def checkSanity(self):
		try:
			if self.__exNum != len(self.__exStarts) or self.__exNum != len(self.__exEnds):
				print(self.__exNum,len(self.__exStarts),len(self.__exEnds))
				raise Exception('the number of exons does not equal the number of exon starting positions or ending position')
			assert(self.__gstart <= self.__gend)
			assert(self.__gstart <= self.__cdstart)
			assert(self.__cdend <= self.__gend)
			for i in range(self.__exNum):
				assert(self.__exStarts[i] <= self.__exEnds[i])
				assert(self.__exStarts[i] >= self.__gstart)
				assert(self.__exEnds[i] <= self.__gend)
		except Exception as err:
			import sys
			print(err)
			sys.exit(0)
	def getRNASeq(self,chr):
		try:
			mrna = ''
			for i in range(self.__exNum):
				if self.__orient=='+':
					mrna = mrna + chr[self.__exStarts[i]:self.__exEnds[i]]
				else:
					mrna = self.reverseCom(chr[self.__exStarts[i]:self.__exEnds[i]]) + mrna		
			return mrna	
		except Exception as err:
			import sys
			print(err)
			sys.exit(0)
	def getCodingRegion(self,chr):
		try:
			cdref = ''
			startp = 0
			for i in range(self.__exNum):
				if startp >= self.__exNum:
					break
				if self.__exStarts[i] <= self.__cdstart and self.__exStarts[startp] < self.__cdstart:
					startp = i
				if self.__exEnds[startp] <= self.__cdstart:
					startp = startp + 1
			while startp < self.__exNum:
				if self.__exStarts[startp] <= self.__cdstart and self.__exEnds[startp] <= self.__cdend:
					if self.__orient == '+':
						cdref = cdref + chr[self.__cdstart:self.__exEnds[startp]]
					else:
						cdref = self.reverseCom(chr[self.__cdstart:self.__exEnds[startp]]) + cdref
				elif self.__exStarts[startp] > self.__cdstart and self.__exEnds[startp] <= self.__cdend:
					if self.__orient == '+':
						cdref = cdref + chr[self.__exStarts[startp]:self.__exEnds[startp]]
					else:
						cdref = self.reverseCom(chr[self.__exStarts[startp]:self.__exEnds[startp]]) + cdref
				elif self.__exStarts[startp] > self.__cdstart and self.__exStarts[startp] < self.__cdend and self.__exEnds[startp] > self.__cdend:
					if self.__orient == '+':
						cdref = cdref + chr[self.__exStarts[startp]:self.__cdend]
					else:
						cdref = self.reverseCom(chr[self.__exStarts[startp]:self.__cdend]) + cdref
				elif self.__exStarts[startp] <= self.__cdstart and self.__exEnds[startp] > self.__cdend:
					if self.__orient == '+':
						cdref = cdref + chr[self.__cdstart:self.__cdend]
					else:
						cdref = self.reverseCom(chr[self.__cdstart:self.__cdend]) + cdref
				startp = startp + 1
			return cdref
			#return startp
			#return self.getRNASeq(chr)	
		except Exception as err:
			import sys
			print(err)
			sys.exit(0)
	def codingRange(self):
		return self.__cdRegion
	def reverseCom(self,seq):
		ref = ''
		for e in seq:
			ref = refGene.dictmap[e] + ref
		return ref	
	def chr(self):
		return self.__chr	
	def orient(self):
		return self.__orient
	def rnaname(self):
		return self.__mrnaid
	def genename(self):
		return self.__alternativeName
	def exlength(self):
		return self.__exlen
	def __str__(self):
		return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t'%(self.__mrnaid,self.__chr,self.__orient,self.__gstart,self.__gend,self.__cdstart,self.__cdend,self.__exNum,','.join(map(str,self.__exStarts))+',',','.join(map(str,self.__exEnds))+',',self.__alternativeName)
	# start from 1	and return [start,end] both end is closed
	def codingRegion(self,quiet = True):
		exons = zip(self.__exStarts,self.__exEnds)
		startEx = 0
		endEx = 0
		for p in range(len(exons)):
			if self.__cdstart >= exons[p][0] and self.__cdstart <= exons[p][1]:
				startEx = p
			if self.__cdend >= exons[p][0] and self.__cdend <= exons[p][1]:
				endEx = p
		
				
		if self.__orient == '+':	
			startpos = sum(self.__exlen[0:startEx]) + self.__cdstart - self.__exStarts[startEx] + 1
			endpos = sum(self.__exlen[0:endEx]) + self.__cdend - self.__exStarts[endEx] 
		else:
			startpos = sum(self.__exlen[endEx + 1:]) -  self.__cdend + self.__exEnds[endEx] + 1
			endpos = sum(self.__exlen[startEx+1:]) + self.__exEnds[startEx] - self.__cdstart   
		if not quiet:
			print(startEx,startpos,endEx,endpos)
		self.__cdRegion = (startpos,endpos)
	def inCodingRegion(self,pos):
		return pos >= self.__cdRegion[0] and pos <= self.__cdRegion[1]	
	# need more work
	def getChrPos(self,mpos):
		if self.__orient == '+':
			off =1 
			index = 0	
			while off <= mpos and index < len(self.__exlen):
				off = off + self.__exlen[index]
				index = index + 1
			if index > 1 and index <= len(self.__exlen):
				extra = mpos -  sum(self.__exlen[0:index-1])
				if extra:
					return self.__exStarts[index -1] + extra 
			else:
				return self.__exStarts[0] + mpos
			
		else:
			off = 1 
			index = len(self.__exlen) - 1
			while off <= mpos and index >= 0:
				off = off + self.__exlen[index]
				index = index - 1
			#extra = mpos - sum(self.__exlen[index+2:])
			extra = off -  mpos
			if index + 1 < len(self.__exlen):
				return self.__exStarts[index+1] +  extra 
	def getMrnaPos(self,pos):
		if self.__exStarts[0] > pos :
			return self.__exStarts[0]
		elif self.__exEnds[-1] < pos:
			return self.__exEnds[-1]
		i = 0
		while True:
			cpos = self.getChrPos(i)
			if cpos == pos:
				return i
			else:
				i = i + 1 
