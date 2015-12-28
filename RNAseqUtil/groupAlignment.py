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
def groupAlignment(mutation,refgeneBymrna,suppfile,outprefix):
	mutationtype = mutation.mutationType()
	start,end = mutation.getChrMutRange()
	norm = []
	supp = []
	if os.path.isfile(suppfile):
		for al in parseSuppFileInd(suppfile):
			mutations = al.getMutations()
			
			mrna = al.getChr()
			for e in  refgeneBymrna[mrna.split('ex')[0]]:
				if e.rnaname() == mrna:
					mrna = e
					break
			checkset = set() 
			for e in mutations:
				m_start = mrna.getChrPos(e[0])
				m_end = mrna.getChrPos(e[1])
				if m_start > m_end:
					m_start,m_end = m_end,m_start
				checkset.add((m_start,m_end,e[2]))	
			if mutationtype == 1 and (start,end,'I') in checkset:
				supp.append(al)
			elif mutationtype == 2 and (start,end,'D') in checkset:
				supp.append(al)
			else:
				norm.append(al)					
	else:
		print('there is no file' + suppfile)
	fhandler = open(outprefix+'.supp','w')
	for e in supp:
		fhandler.write(e.__str__())
	fhandler.close()
	fhandler = open(outprefix+'.norm','w')
	for e in norm:
		fhandler.write(e.__str__())
	fhandler.close()
		
# given a set of suppfiles and a mutation set
# return all supp alignment that cover the mutations in the mutation set 
def extractMappedImp(suppfiledir,mutations,refgeneBymrna,outprefix):
		
	for m in mutations:
		chr = m.getChromo()
		start,end=m.getChrMutRange()
		filename = '%s_%s-%s'%(chr,start,end)
		absfilename = os.path.join(suppfiledir,filename)+'.supp'
		tmpprefix = os.path.join(outprefix,filename)	
		groupAlignment(m,refgeneBymrna,absfilename,tmpprefix)	
			
def extractMapped(sample,mutations,refGeneByMrna,outdir):
	suppdir = '/home/luzhao/THR104_analysis/all_alignment'
	suppfiledir = os.path.join(suppdir,sample)
	outprefix =os.path.join(outdir,sample)
	if not os.path.isdir(outprefix):
		os.mkdir(outprefix)
	print(sample,suppfiledir,outprefix)
	extractMappedImp(suppfiledir,mutations,refGeneByMrna,outprefix)
if __name__ == '__main__':
	from multiprocessing import Process
	import time
	outdir = '/home/luzhao/THR104_analysis/suppAlignment'
	mutPrefix = '/home/luzhao/THR104_analysis/finalmutation/mutation_common_except_normal'
	data_infor = sampleInfo()
	samples = data_infor.samples()
	print(samples)
	samplelist = []
	for k,v in samples.items():
		samplelist += v

	sys.stdout.flush()
	refgenebymrna = parseRefGene('/home/luzhao/THR104_analysis/chromo/refGene.txt')
	mutations = parseMutation(mutPrefix + '.csv',mutPrefix+'.supp',mutPrefix + '.norm')
	proc = []
	for s in samplelist:
		p = Process(target = extractMapped,args=(s,mutations,refgenebymrna,outdir))
		proc.append(p)
		p.start()
	for e in proc:
		e.join()
	print('finish the extracting the all alignments')	
	
