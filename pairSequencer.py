#!/usr/bin/python

from sequenceError import sequenceError
from reads import read
from sam import Sam
from fragment import Fragment
from singleSequencer import singleSequencer

#given a fragment,generate a read from it

# mimic the sequencing process that happen at a single fragment for the single end reads
# The process are described as follows:
# 1. the each read is generated from the prime position
# 2.currently, the reads is generated from 
# 3. repeat step 1
class pairSequencer(singleSequencer):
	def __init__(self,errorModel,rdlen):
		singleSequencer.__init__(self,errorModel,rdlen)
	
	# given the mrna sequence, the start and expand of the fragment at this mrna
	def generate(self,frag,sequence,rdID):
		length = len(frag)
		start = frag.startPos()
		if length < 2*self.rlen:
			return (None,None) 
		else:
			res = []
			rd,align = self.generateImp(start,sequence,frag.transcriptID(),rdID+'_1')
			res.append((rd.__str__(),align.__str__()))
			start = start + length - self.rlen
			rd,align = self.generateImp(start,sequence,frag.transcriptID(),rdID+'_2',True)
			rd.complem()
			align.setFlag(16)
			res.append((rd.__str__(),align.__str__()))			
			return res

if __name__ == '__main__':
	sequence = 'ACGACGACGAGCCCGGGTTTTCAGG'
	model = sequenceError()
	sequencing = pairSequencer(model,5)
	frag = Fragment('hello',0,len(sequence))
	print(sequencing.generate(frag,sequence,'simulated'))
					
