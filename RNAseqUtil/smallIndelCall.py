#!/usr/bin/python

from multiprocessing import Queue,Process
from validateMutation import checkAlign
from mutation import mutation
from mutation_detect import loadPartialSamByPos
from mutation_detect import getMappedCoverage
from alignment import alignment
import sys
import os
import threading,time
class smallIndelCall(Process):
	def __init__(self,samfile,outdir,refmrna,refgenes,cutoff,inputmrna,prefix):
		self.samfile = samfile
		self.outdir = outdir
		self.prefix = prefix
		self.refmrna = refmrna
		self.refgenes = refgenes
		self.cutoff = cutoff
		self.inputmrna = inputmrna
		Process.__init__(self)
	def run(self):
		print('start the small indel procedure')
		sys.stdout.flush()
	#	dictsam = loadPartialSamByPos(self.samfile,self.inputmrna)
#		for i in self.inputmrna: 
		mutations = self.mutationSam(self.samfile,self.inputmrna)
		validMut = self.filterMutations(mutations)
		mut,supp,norm = self.outputMutations(validMut)
		if mut and supp:
			print(self.prefix)
			mutfile = os.path.join(self.outdir,self.prefix + '.csv')
			suppfile = os.path.join(self.outdir,self.prefix + '.supp')
			normfile = os.path.join(self.outdir,self.prefix + '.norm')
			print(mutfile,suppfile,normfile)
			muth = open(mutfile,'w')
			supph = open(suppfile,'w')
			normh = open(normfile,'w')
			muth.write(mut)
			supph.write(supp)
			normh.write(norm)
			muth.close()
			supph.close()
			normh.close()	
		print('process %s finished' % (os.getpid()))
	def mutationSam(self,samfile,inputmrna):
        	try:
#			import time
               		mutations  = {}
			handler = open(self.samfile,'r')
			total = 0
#			queritime = 0
#			totaltime = 0
#			checkAlignTime = 0
#			updatetime = 0
			for line in handler:
#				totalc = time.time()
				if line[0] == '@':
					continue
				total = total + 1
				rec = line.rstrip('\n').split()
#				c = time.time()
				if not rec[2] in inputmrna:
					continue
				#lc = time.time()
				#queritime += (lc - c)
				if not total %1000000:
					print('%s finished %s\t%s ' %(os.getpid(),total,len(mutations)))
					sys.stdout.flush()
				align = alignment(rec[0],rec[2],rec[1],rec[3],rec[9],rec[5])
				#totala = time.time()	
                               	muts = checkAlign(align,self.refmrna)
				#totalal = time.time()
				#checkAlignTime += (totalal - totala)
                               	for e in muts:
                                 	if e.__hash__() in mutations:
                                             	mutations[e.__hash__()].addSupport(e.getSupport())
                                       	else:
                                               	mutations[e.__hash__()] = e
				#updatetime += (time.time() - totalal)
				#totall = time.time()
				#totaltime += (totall - totalc)
                	return mutations
        	except Exception as err:
                	print(err)
                	sys.exit(0)
	def filterMutations(self,mutations):
        	validmut = {}
		dictsam = {}
        	if mutations:
                	for k,v in mutations.items():
                               	rnas = self.refgenes[v.mutationChr().split('ex')[0]]
                               	for e in rnas:
                                       	if e.rnaname() == v.mutationChr():
                                               	if e.inCodingRegion(v.mutationPos()):
                                                       	v.setChromosome(e.chr())
                                                       	v.setChrMutRange(e)
                                                       	v.setGene(e.genename())
                                                       	normsupport = getMappedCoverage(v,dictsam)
                                                       	v.addNormalSupport(normsupport)
							key = str(v.mutationType()) + v.getChromo() + ':' + str(v.getChrMutRange()[0]) + v.mutationSeq()
							if key in validmut:
                                                       		validmut[key].addSupport(v.getSupport())
							else:
								validmut[key] = v
		res = []
		for k,v in validmut.items():
			if v.NumOfSupportReads() >= self.cutoff:
				res.append(v)
        	return res
	def outputMutations(self,mutations):
        	mutfile = ''
        	suppfile = ''
        	normfile = ''
        	if mutations:
                	mutations = list(set(mutations))
                	mutations.sort(reverse=True)
                	for e in mutations:
                        	mutfile = mutfile + e.__str__()+'\n'
                        	for s in e.getSupport():
                                	suppfile = suppfile + s.__rec__()+'\n'
                        	for s in e.getNormalSupport():
                                	normfile = normfile + s.__rec__()+'\n'
		return (mutfile,suppfile,normfile)
