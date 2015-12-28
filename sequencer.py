#!/usr/bin/python

from sequenceError import sequenceError
from reads import read
from sam import Sam
from fragment import Fragment
from singleSequencer import singleSequencer
from pairSequencer import pairSequencer

#given a fragment,generate a read from it

# mimic the sequencing process that happen at a single fragment for the single end reads
# The process are described as follows:
# 1. the each read is generated from the prime position
# 2.currently, the reads is generated from 
# 3. repeat step 1
class Sequencer(object):
	def __init__(self,errorModel,rdlen,pair = False):
		if pair:
			self.sequencer = pairSequencer(errorModel,rdlen)
		else:
			self.sequencer = singleSequencer(errorModel,rdlen)
	# given the mrna sequence, the start and expand of the fragment at this mrna
	def generate(self,frag,sequence,rdID):
		return self.sequencer.generate(frag,sequence,rdID)

if __name__ == '__main__':
	sequence = 'ACGACGACGAGCCCGGGTTTTCAGG'
	model = sequenceError()
	sequencing = Sequencer(model,5)
	frag = Fragment('hello',0,10)
	print(sequencing.generate(frag,sequence,'simulated'))
					
