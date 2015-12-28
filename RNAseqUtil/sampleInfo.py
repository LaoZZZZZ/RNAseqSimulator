#!/usr/bin/python
import glob
import sys
import os
class sampleInfo:
	directory = '/home/stc/Documents/GenomicData/RNAseq/WBahou'	
	sampleReads = {'A019':'multiple','A023':'multiple','N084':'multiple','THR039':'multiple','THR053':'single','THR100':'single','THR101':'single','THR102':'single','THR104':'single','THR116':'single','THR136':'single','THR137':'single'}
	def samples(self):
		return {'normal':['A019','A023','N084','THR039','THR053'],'ET':['THR100','THR101','THR102','THR104','THR116','THR136','THR137']}
	def readsData(self,sample):
		os.chdir(os.path.join(self.directory, sample))
		if sampleInfo.sampleReads[sample] == 'multiple':
			return glob.glob('*[0-9].gz')
		else:
			return glob.glob('*sequence.gz')	
	def ET(self):
		return ['THR100','THR101','THR102','THR104','THR116','THR136','THR137']
	def norm(self):
		return 	['A019','A023','N084','THR039','THR053']		
if __name__ == '__main__':
	samples = sampleInfo()
	print(samples.readsData('A019'))		
	print(samples.readsData('THR102'))
	print(samples.readsData('THR039'))
	
