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
# given a mrna object and a region on this mrna
# return the corresponding range on the genome reference
def getGenomePos(refgene,start,end):
	chrStart = refgene.getChrPos(start)
	chrEnd = refgene.getChrPos(end)
	if chrStart <= chrEnd:
		return (chrStart,chrEnd)
	else:
		return (chrEnd,chrStart)
def getGenomeCoverage(align,refgenes)
	start,end = align.getRange()
	mrna = align.getChr()
	refgene = refgenes[mrna.split('ex')[0]]
	for e in refgene:
		if e.rnaname() == mrna:
			mrna = e
			break
	return getGenomePos(mrna,start,end) 

def extractMappedImp(samfile,mutations,refgene,outprefix):
	rnas = {}
	fhandlers = {}
	coverRanges = set()
	for m in mutations:
		for rna in refgene[m.getGene()]:
			rnas[rna.rnaname()] = []
			filename = os.path.join(outprefix,rna.rnaname()+'.supp')
			handler = open(filename,'w')
			fhandlers[rna.rnaname()] = handler
			
	samhandler = SESamfileIterator(samfile)
	 
	for line in samhandler:
		if line.getChr() in rnas:
			rnas[line.getChr()].append(line)
			if len(rnas[line.getChr()]) > 1000:
				for e in rnas[line.getChr()]:
					fhandlers[line.getChr()].write(line.__str__())
				rnas[line.getChr()] = []
	
	for k,v in rnas.items():
		for e in v:
			fhandlers[k].write(e.__str__())
	for k,v in fhandlers:
		v.close()			
def extractMapped(sample,mutPrefix,refGene,outdir):
	samdir = '/home/luzhao/THR104_analysis/mapped'
	mutationdir = '/home/luzhao/THR104_analysis/finalmutation'
	samfile = os.path.join(samdir,sample+'.sam')
	mutations = parseMutation(os.path.join(mutationdir,mutPrefix + '.csv'),os.path.join(mutationdir,mutPrefix+'.supp'),os.path.join(mutationdir,mutPrefix + '.norm'))
	outprefix =os.path.join(outdir,sample)
	if not os.path.isdir(outprefix):
		os.mkdir(outprefix)
	extractMappedImp(samfile,mutations,refGene,outprefix)
if __name__ == '__main__':
	from multiprocessing import Process
	import time
	outdir = '/home/luzhao/THR104_analysis/mutationSam'
	data_infor = sampleInfo()
	samples = data_infor.samples()
	print(samples)
	samplelist = []
	for k,v in samples.items():
		samplelist += v

	sys.stdout.flush()
	refgene = parseRefGeneByGene('/home/luzhao/THR104_analysis/chromo/refGene.txt')
	proc = []
	sz = 0
	for s in samplelist:
		p = Process(target = extractMapped,args=(s,'mutation_common_except_normal',refgene,outdir))
		proc.append(p)
		p.start()
		if sz >= 2:
			for e in proc:
				e.join()
			proc = []
			sz = 0		
		sz += 1
	for e in proc:
		e.join()
	print('finish the extracting the all alignments')	
	
