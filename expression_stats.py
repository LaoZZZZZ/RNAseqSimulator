#!/usr/bin/python


from RNAseqUtil.samIterator import SESamfileIterator as samIter
from RNAseqUtil.alignment import alignment
from mrnaGenerator import mrnaGenerator
# calculate the expression level of each transcript using 
class expression_stats(object):
	def __init__(self,expression,samfile,outfile):
		self.outfile = outfile
		self.express = {}
		for k,v in expression.items():
			self.express[k] = list(v) + [0] 	
		self.parseSam(samfile)
		self.stats()
	def parseSam(self,samfile):
		self.total = 0.0
		self.rdlen = 0.0
		print(samfile)
		for al in samIter(samfile):
			if al:
				self.rdlen = len(al)
				break
		for al in samIter(samfile):
			if al:
				self.express[al.getChr()][-1] += 1
				self.total += 1
	def stats(self):
		self.normalized_reads = 0.0
		if self.total:
			for k,v in self.express.items():
				rpkm = v[-1] * 10**9/(v[0] * self.total)
				tpm = v[-1]*10**6/v[0] 
				self.normalized_reads += float(v[-1])/v[0]
				v += [rpkm,tpm]
			for k,v in self.express.items():
				v[-1] = v[-1] / self.normalized_reads
	def dump(self):
		header = 'transcriptID,length,absolute abundance, relative abundance,expected read count,RPKM,TPM\n'
		handler = open(self.outfile,'w')
		handler.write(header)
		for k,v in self.express.items():
			handler.write('%s,%s\n'%(k,','.join(map(str,v))))
		handler.close()
					
