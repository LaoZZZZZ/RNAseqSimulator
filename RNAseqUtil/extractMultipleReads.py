#!/usr/bin/python
import os
import sys
from extract_mapped import extractRegionMapping
from parseRefGene import parseRefGeneByGene
from refGeneStruct import refGene	
from samIterator import SESamfileIterator
from parseMutSupp import parseMutation
from mutation import mutation
from alignment import alignment
from sampleInfo import sampleInfo
def extractMultipleImp(samfile,mutations,outprefix):
	readsIDs = set() 
	aligns = []
	for m in mutations:
		matches = m.getSupport()
		for a in matches:
			id = a.getID()
			readsIDs.add(id)
		norm = m.getNormalSupport()
		for a in norm:
			id = a.getID()
			readsIDs.add(id)
	samhandler = SESamfileIterator(samfile)
	for line in samhandler:
		if line.getID() in readsIDs:
			aligns.append(line.__str__())
	fhandler = open(outprefix+'_mutiple.csv','w')	
	for k in sorted(aligns):
		fhandler.write(k)
	fhandler.close()			
def extractMultipleReads(sample,mutPrefix,outdir):
	samdir = '/home/luzhao/THR104_analysis/mapped'
	mutationdir = '/home/luzhao/THR104_analysis/finalmutation/possible_mutations'
	samfile = os.path.join(samdir,sample+'.sam')
	mutations = parseMutation(os.path.join(mutationdir,mutPrefix + '.csv'),os.path.join(mutationdir,mutPrefix+'.supp'),os.path.join(mutationdir,mutPrefix + '.norm'))
	outprefix = os.path.join(outdir,sample)
	extractMultipleImp(samfile,mutations,outprefix)
if __name__ == '__main__':
	from multiprocessing import Process
	import time
	outdir = '/home/luzhao/THR104_analysis/multipleAlignment'
	data_infor = sampleInfo()
	samples = data_infor.samples()
	print(samples)
	samplelist = []
	for k,v in samples.items():
		samplelist += v

	sys.stdout.flush()
	#extractMultipleReads(samples['ET'][0],samples['ET'][0]+'_small',outdir)
	proc = []
	for s in samplelist:
		p = Process(target = extractMultipleReads,args=(s,s+'_small',outdir))
		proc.append(p)
		p.start()
	for e in proc:
		e.join()
	print('finish the extracting the all alignments')	
	
