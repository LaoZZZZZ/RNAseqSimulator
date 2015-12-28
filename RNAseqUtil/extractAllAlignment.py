#!/usr/bin/python
import os
import sys
from extract_mapped import extractRegionMapping
from parseRefGene import parseRefGeneByGene,parseRefGene
from refGeneStruct import refGene	
from samIterator import SESamfileIterator
from parseMutSupp import parseMutation,parseSuppFileInd
from mutation import mutation
from alignment import alignment
from sampleInfo import sampleInfo
from mrna2genome import checkCoverage
# given a mutation range, a supp file, and refgenes db
# extract all alignments in suppfile that cover that mutation range
# return a list of these alignments
def findOverlapAlignment(mutRange,suppfile,refgenes):
	chr,mut = mutRange.split(':')
	mut_start,mut_end = mut.split('-')
	mut_start = int(mut_start)
	mut_end = int(mut_end)
	result = []
	if os.path.isfile(suppfile):
		for al in parseSuppFileInd(suppfile):
			if checkCoverage(al,refgenes,chr,mut_start,mut_end):
				result.append(al)
	else:
		print('there is no file' + suppfile)
	return result	
		
# given a set of suppfiles and a mutation set
# return all supp alignment that cover the mutations in the mutation set 
def extractMappedImp(suppfiledir,mutations,refgeneBygene,refgeneBymrna,outprefix):
	coverRanges = {}
	for m in mutations:
		chr = m.getChromo()
		start,end=m.getChrMutRange()
		key = '%s:%s-%s'%(chr,start,end)
		coverRanges[key] = []
		for m in  refgeneBygene[m.getGene()]:
			coverRanges[key].append(m.rnaname())
	for k,v in coverRanges.items():
		aligns = []
		chr,range = k.split(':')
		name = '%s_%s'%(chr,range)
		fhandler = open(os.path.join(outprefix,name+'.supp'),'w')
		for f in v:
			aligns += findOverlapAlignment(k,os.path.join(suppfiledir,f+'.supp'),refgeneBymrna)	
		print(len(aligns))
		for a in aligns:
			fhandler.write(a.__str__())
		fhandler.close()
			
def extractMapped(sample,mutations,refGeneByGene,refGeneByMrna,outdir):
	suppdir = '/home/luzhao/THR104_analysis/mutationSam'
	suppfiledir = os.path.join(suppdir,sample)
	outprefix =os.path.join(outdir,sample)
	if not os.path.isdir(outprefix):
		os.mkdir(outprefix)
	print(sample,suppfiledir,outprefix)
	extractMappedImp(suppfiledir,mutations,refGeneByGene,refGeneByMrna,outprefix)
if __name__ == '__main__':
	from multiprocessing import Process
	import time
	outdir = '/home/luzhao/THR104_analysis/all_alignment'
	mutPrefix = '/home/luzhao/THR104_analysis/finalmutation/mutation_common_except_normal'
	data_infor = sampleInfo()
	samples = data_infor.samples()
	print(samples)
	samplelist = []
	for k,v in samples.items():
		samplelist += v

	sys.stdout.flush()
	refgenebygene = parseRefGeneByGene('/home/luzhao/THR104_analysis/chromo/refGene.txt')
	refgenebymrna = parseRefGene('/home/luzhao/THR104_analysis/chromo/refGene.txt')
	mutations = parseMutation(mutPrefix + '.csv',mutPrefix+'.supp',mutPrefix + '.norm')
	proc = []
	for s in samplelist:
		p = Process(target = extractMapped,args=(s,mutations,refgenebygene,refgenebymrna,outdir))
		proc.append(p)
		p.start()
	for e in proc:
		e.join()
	print('finish the extracting the all alignments')	
	
