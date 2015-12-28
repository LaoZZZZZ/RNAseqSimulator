#!/usr/bin/python

from sequenceError import sequenceError
from reads import read
from sam import Sam
from fragment import Fragment


#given a fragment,generate a read from it

# mimic the sequencing process that happen at a single fragment for the single end reads
# The process are described as follows:
# 1. the each read is generated from the prime position
# 2.currently, the reads is generated from 
# 3. repeat step 1
class singleSequencer(object):
	def __init__(self,errorModel,rdlen):
		self.errorModel = errorModel
		self.rlen = rdlen
	# given the mrna sequence, the start and expand of the fragment at this mrna
	def generate(self,frag,sequence,rdID):
		length = len(frag)
		start = frag.startPos()
		if length < self.rlen:
			return (None,None) 
		else:
			rd,align = self.generateImp(start,sequence,frag.transcriptID(),rdID)
			return [(rd.__str__(),align.__str__())]
	# the actual sequence action
	def generateImp(self,start,sequence,seqID,rdID,rev=False):
		self.cigar = ''		
		seq = ''
		qua = ''
		points = range(self.rlen)
		if rev:
			points.reverse()
		for off in range(len(points)):
			b,q= self.errorModel.basePair(sequence[start+points[off]],off)	 
			seq += b
			qua += chr(q)
			if b != sequence[start+off]:
				self.cigar += 'X'
			else:
				self.cigar += 'M'
		rd = read(rdID,seq,qua)
		alignment = Sam(rdID,start+1,seqID,self.generateCigar())
		alignment.setSeq(seq)
		alignment.setQua(qua)
		alignment.setMapQ(0)
		return (rd,alignment)
	def generateCigar(self):
		cigar = ''
		tmp = ''	
		for ch in self.cigar:
			if not tmp or ch == tmp[-1]:
				tmp += ch
			else:
				cigar += str(len(tmp))
				cigar += tmp[0]
				tmp = ch	
		if tmp:
			cigar += str(len(tmp))
			cigar += tmp[0]
		return cigar			

if __name__ == '__main__':
	sequence = 'ACGACGACGAGCCCGGGTTTTCAGG'
	model = sequenceError()
	sequencing = singleSequencer(model,5)
	frag = Fragment('hello',0,10)
	print(sequencing.generate(frag,sequence,'simulated'))
					
