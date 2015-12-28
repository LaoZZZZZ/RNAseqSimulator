#!/usr/bin/python
from multiprocessing import Process
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
from mutationSet import mutationSet
import glob
import shutil
#list all mutations that happened in normal sample
def falsePositive(norm_samples,srcdir,outdir):
	falsePositive =  set()
	for s in norm_samples:
		directory = os.path.join(srcdir,s)
		os.chdir(directory)
		if not os.path.isdir(os.path.join(outdir,s)):
			os.makedirs(os.path.join(outdir,s))
		for e in glob.glob('*.supp'):
			if os.path.getsize(e) > 0:
				shutil.copy(e,os.path.join(outdir,s)+'/' + e)
				shutil.copy(e.split('.')[0]+'.norm',os.path.join(outdir,s)+'/' + e.split('.')[0]+'.norm')
				falsePositive.add(e.split('.')[0])
	return falsePositive 
def updateSingle(mutations,sample,outdir,fp,possibledir):
	mutations = mutationSet([os.path.join(inputmutdir,sample + '_small.csv'),os.path.join(inputmutdir,sample + '_small.supp'),os.path.join(inputmutdir,sample + '_small.norm')])
	if not os.path.isdir(os.path.join(outdir,sample)):
		os.makedirs(os.path.join(outdir,sample))
	if not os.path.isdir(os.path.join(possibledir,sample)):
		os.makedirs(os.path.join(possibledir,sample))
	muts = mutations.getMutations()
	mutset = {}
	possmutset = {}
	for k,v in muts.items():
		chr = v.getChromo()
		start,end = v.getChrMutRange()
		filename = '%s_%s-%s'%(chr,start,end)
		if filename in fp:
			possmutset[k] = v
		else:
			mutset[k] = v
	possiblePrefix = os.path.join(possibledir,sample+'_small')
	mut,supp,norm = mutationSet(possmutset).__str__()
	open(possiblePrefix+'.csv','w').write(mut)
	open(possiblePrefix+'.norm','w').write(norm)
	open(possiblePrefix+'.supp','w').write(supp)
	mutPrefix = os.path.join(outdir,sample+'_small')
	mut,supp,norm = mutationSet(mutset).__str__()
	open(mutPrefix+'.csv','w').write(mut)
	open(mutPrefix+'.norm','w').write(norm)
	open(mutPrefix+'.supp','w').write(supp)
		
	
			
# given a set of suppfiles and a mutation set
# return update the support alignment and normal alignment 
def update(suppfiledir,inputmutdir,norm_samples,ET_samples,possibleMutoutdir,outdir):
	fp = falsePositive(norm_samples,suppfiledir,possibleMutoutdir)
	proc = []
	for s in ET_samples:
		print(s)
		p = Process(target = updateSingle,args = (inputmutdir,s,outdir,fp,possibleMutoutdir))
		p.start()
		proc.append(p)
	for e in proc:
		e.join()
def loadMutationsByLine(file):
	mutations = {}
	fhandler = open(file,'r')
	for line in fhandler:
		if line[0]=='@':
			continue
		else:
			rec = line.rstrip().split(',')
			key = rec[3]+'_'+rec[4]+'-'+rec[5]
			mutations[key] = rec[0:15]
	return mutations
def addAlignment(sample,mutations,suppfiledir,outdir):
	os.chdir(os.path.join(suppfiledir,sample))
	final = {}
	suppal = ''
	normal = ''
	for k,v in mutations.items():
		suppfile = k+'.supp'
		normfile = k+'.norm'
		supp = set()
		numSupp = 0
		suppID = ''	
		for al in parseSuppFileInd(suppfile):
			if not al.getSeq() in supp:	
				supp.add(al.getSeq())
				numSupp += 1
			if suppID != '':
				suppID += (' ' + al.getID())
			else:
				suppID = al.getID()
		print(numSupp)
		suppID.rstrip()
		norm = set()
		numNorm = 0
		normID = ''
		for al in parseSuppFileInd(normfile):
			if not al.getSeq() in norm:
				numNorm += 1
				norm.add(al.getSeq())
			if normID:
				normID += (' ' + al.getID())
			else:
				normID = al.getID()
		mut = v
		mut[9] = numSupp
		mut[10] = len(suppID.split())
		mut[11] = numNorm
		mut[12] = len(normID.split())
		if numSupp > 0:
			mut[13] = float(numSupp)/(numSupp + numNorm)
			mut[14] = float(mut[10])/(mut[10] + mut[12])
		else:
			mut[13] = 0
			mut[14] = 0
		mut[9:15] = list(map(str,mut[9:15]))
		mut.append(suppID)
		mut.append(normID)
		final[k] = ','.join(mut)
		suppal += open(suppfile).read()
		normal += open(normfile).read()
	fhandler = open(os.path.join(outdir,sample+'_small.csv'),'w')
	for k, v in final.items():
		fhandler.write(v+'\n')
	fhandler.close()
	fhandler = open(os.path.join(outdir,sample+'_small.supp'),'w')
	fhandler.write(suppal)
	fhandler.close()
	fhandler = open(os.path.join(outdir,sample+'_small.norm'),'w')
	fhandler.write(normal)
	fhandler.close()
									
def extractMapped(sample,mutations,refGeneByMrna,outdir):
	suppdir = '/home/luzhao/THR104_analysis/all_alignment'
	suppfiledir = os.path.join(suppdir,sample)
	outprefix =os.path.join(outdir,sample)
	if not os.path.isdir(outprefix):
		os.mkdir(outprefix)
	print(sample,suppfiledir,outprefix)
	extractMappedImp(suppfiledir,mutations,refGeneByMrna,outprefix)
if __name__ == '__main__':
	import time
	suppfiledir = '/home/luzhao/THR104_analysis/suppAlignment'
	inputmutdir = '/home/luzhao/THR104_analysis/finalmutation/updated'
	mutPrefix = '/home/luzhao/THR104_analysis/finalmutation/'
	possibleMutoutdir = '/home/luzhao/THR104_analysis/possibleMutation'
	outdir = '/home/luzhao/THR104_analysis/finalmutation/final'
	data_infor = sampleInfo()
	samples = data_infor.samples()
	print(samples)
	samplelist = []
	for k,v in samples.items():
		samplelist += v

	sys.stdout.flush()
	#update(suppfiledir,inputmutdir,samples['normal'],samples['ET'],possibleMutoutdir,outdir)
	proc = []
	mutations = loadMutationsByLine(os.path.join(outdir,'THR100_small.csv'))
	for s in samples['normal']:
		print(s)
		p = Process(target = addAlignment,args=(s,mutations,suppfiledir,outdir))
		proc.append(p)
		p.start()
	for e in proc:
		e.join()
	print('finish the extracting the all alignments')	
	
