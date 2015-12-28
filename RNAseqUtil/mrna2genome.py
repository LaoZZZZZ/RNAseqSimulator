#!/usr/bin/python
import os
import sys
from parseRefGene import parseRefGene
from refGeneStruct import *	
from mutation import mutation
from alignment import alignment
from sampleInfo import sampleInfo
from parseMutSupp import parseSuppFileInd
def getGene(align,refgenes):
	mrna = align.getChr()
	if mrna.count('_') >= 2:
		rec = mrna.split('_')
		mrna = '%s_%s'%(rec[0],rec[1])
	key =mrna.split('ex')[0]
	if not key in refgenes:
		return None
	refgene = refgenes[key]
	for e in refgene:
		if e.rnaname() == mrna:
			mrna = e
			break
	if isinstance(mrna,refGene):
		return mrna.genename()
	else:
		return None
		raise ValueError('can not find genes for %s'%(mrna))
def getAlignmentChr(align,refgenes):
	mrna = align.getChr()
	refgene = refgenes[mrna.split('ex')[0]]
	for e in refgene:
		if e.rnaname() == mrna:
			mrna = e
			break
	return mrna.chr()
# given a mrna object and a region on this mrna
# return the corresponding range on the genome reference
def getGenomePos(refgene,start,end):
	chrStart = refgene.getChrPos(start)
	chrEnd = refgene.getChrPos(end)
	if chrStart <= chrEnd:
		return (chrStart,chrEnd)
	else:
		return (chrEnd,chrStart)
# given a alignment on a mrna
# return the alignmetn coverage on genome
def getGenomeCoverage(align,refgenes):
	mrna = align.getChr()
	refgene = refgenes[mrna.split('ex')[0]]
	for e in refgene:
		if e.rnaname() == mrna:
			mrna = e
			break
	start,end = align.getRange()
	return getGenomePos(mrna,start,end) 
# check two range is overlapped
def isOverlaped(start1,end1,start2,end2):
	if end1 < start2 or start1 > end2:
		return False
	return True
# check if an alignment cover the range of the mutation range
# not necessary the read carries the mutation signature
def checkCoverage(align,refgenes,chr,mut_start,mut_end):
	if(getAlignmentChr(align,refgenes) == chr):
		al_start,al_end = getGenomeCoverage(align,refgenes)
		return isOverlaped(al_start,al_end,mut_start,mut_end)
	return False
if __name__ == '__main__':
	direct = '/home/luzhao/THR104_analysis/mutationSam/A019'
	refgenes = parseRefGene('/home/luzhao/THR104_analysis/chromo/refGene.txt') 
	import os
	import sys
	os.chdir(direct)
	import glob
	for f in glob.glob('*.supp'):
		for e in parseSuppFileInd(f):
			print(e)
			print(getGenomeCoverage(e,refgenes))
			sys.exit()	
