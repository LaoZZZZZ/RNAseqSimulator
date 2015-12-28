#!/usr/bin/python
import os
import sys
from refGeneStruct import refGene
from parseRefGene import parseRefGene
from mutation_detect import loadPartialSamByPos
from smallIndelCall import smallIndelCall
from multiprocessing import Process,Queue
import random
def splitMrna(mrna,n):
	size = len(mrna)/n
	left = len(mrna)%n
	split = []
	random.shuffle(mrna)
	for i in range(n):
		sub = mrna[i*size:(i+1)*size]
		split.append(sub)
	if left:
		split[-1] = split[-1] + mrna[n*size:]
	return split
def samDrive(samfile,genomefile,origsamfile,outdir,prefix,cutoff,refgfile):
	refgenes = parseRefGene(refgfile)
	mutd = mutationSam(samfile,genomefile)
	mutations = filterMutations(mutd,cutoff,refgenes,origsamfile)
	outputMutations(outdir,prefix,mutations)	
def samDriveOnline(samfile,refmrna,outdir,prefix,cutoff,refgenes,nproc,split):
	import time
	print('start sam drive')
	sys.stdout.flush()
	
	p = smallIndelCall(samfile,outdir,refmrna,refgenes,cutoff,split,prefix)
	p.start()
	p.join()
	print('finished sam drive')
				
if __name__ == '__main__':
	dir = '/home/luzhao/THR104_analysis'
	refgenes = parseRefGene(os.path.join(dir,'chromo/refGene.txt'))
	#inputdir = os.path.join(dir,'mutation')
	outdir = os.path.join(dir,'finalmutation')
	origsamfile = os.path.join(dir,'THR104.sam')
	#extractCodingMutDrive(refgene,inputdir,outdir,origsamfile)
	refmrna = loadRefMrna(os.path.join(dir,'chromo/ref/mrna.fa'))
	inputdir = os.path.join(dir,'mapped')
	samfile = os.path.join(dir,'mapped/test.sam')

	samDriveOnline(samfile,refmrna,origsamfile,outdir,'mrna',1,refgenes)
#	samDriveOnline(samfile,ref,origsamfile,os.path.join(dir,'mapped'),'mrna',2,refgene)
	#mutationSam(samfile,ref)
