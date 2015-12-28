#!/usr/bin/python

from fragmentGenerator import fragmentGenerator as fGenerator

from mrnaDetermineGenerator import mrnaDetermineGenerator
from amplification import amplifier
from reads import read
from mrnalibrary import mrnaLibrary
from sequenceError import sequenceError
from sequencer import Sequencer
from sequenceError import sequenceError
# given the mrna template library and fragment generator
# generat reads and output them into the specifed file and format
# if necessary, the alignment result will be output to the specified file in samfile format
class readGenerator(object):
	def __init__(self,machine,mrnalib,fgenerator,outprefix,alignment=False):
		self.mrnalib = mrnalib
		self.fgen = fgenerator
		self.outprefix = outprefix
		self.outAlignment = alignment
		self.seqMachine = machine
		self.size= 0
	def simulatorDrive(self):
		print('In the process of sequencing')
		self.buffer = []
		total = 0
		readsFile = open(self.outprefix+'.fq','w')	
		samfile = open(self.outprefix+'.sam','w')
		for f in self.fgen:
			#print(f)
			chr = f.transcriptID()
			seq = self.mrnalib[chr]
			res = self.seqMachine.generate(f,seq[1],'simulated'+str(total))
			if res:
				for rd,sam in res:	
					readsFile.write(rd+'\n')
					samfile.write(sam+'\n')
				total += 1
		print('Sequence finished')
		readsFile.close()
		samfile.close()
		self.size = total
	def __del__(self):
		print('totally generated %s reads'%(self.size))		



if __name__ == '__main__':
	mrnafile = 'mrna.txt'  #template library
	rdlen = 50      #read length
	poolsize = 10  #size of transcriptome set
	lambd = 120  # fragment rate
	length_lower = 100 
	length_upper = 140 

	mrnalib = mrnaLibrary(mrnafile) #object that holds all template 

	mrnaGen = mrnaDetermineGenerator(mrnafile,poolsize)	
	fragGen = fGenerator(mrnaGen,lambd,length_lower,length_upper)  # generate all fragments given the transcriptome
	
	outpref = 'RNASimulation'  # out put prefix
	ErrorModel = sequenceError()  # define the  sequence error model
	seqmachine = Sequencer(ErrorModel,rdlen)
	
	simulator = readGenerator(seqmachine,mrnalib,fragGen,outpref)
	simulator.simulatorDrive() 
			 
