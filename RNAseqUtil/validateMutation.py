#!/usr/bin/python

from refGeneStruct import refGene
from parseRefGene import parseRefGene
from parseMutSupp import parseMutation
from mutation import mutation
from align_candidate import align_candidate
from alignment import alignment
from strMatch import *
from reads import read
from mutation_detect import getMappedCoverage
from mutation_detect import loadSamByPos
import os
import sys
import re
from parseRefGene import loadRefMrna
from parseRefGene import parseRefGene
from extract_mapped import extractSamRegion
def outputMutations(outdir,prefix,mutations):
	mutfile = os.path.join(outdir,prefix+'.csv')
	suppfile = os.path.join(outdir,prefix+'.supp')
	normfile = os.path.join(outdir,prefix+'.norm')	
	print(mutfile,suppfile)
	if mutations:
		muth = open(mutfile,'w')
		supph = open(suppfile,'w')
		normh = open(normfile,'w')
		mutations = list(set(mutations))
		mutations.sort(reverse=True)
		for e in mutations:
			muth.write(e.__str__()+'\n')
			for s in e.getSupport():
				supph.write(s.__rec__()+'\n')
			for s in e.getNormalSupport():
				normh.write(s.__rec__()+'\n')	
		muth.close()
		supph.close()
		normh.close()
	
def filterMutations(mutations,cutoff,refgene,dictsam):
	validmut = []
	if mutations:
		for k,v in mutations.items():
			rnas = refgene[v.mutationChr().split('ex')[0]]
			for e in rnas:
				if e.rnaname() == v.mutationChr():
					if e.inCodingRegion(v.mutationPos()):
						v.setChromosome(e.chr())
						v.setChrMutRange(e)
						v.setGene(e.genename())
						normsupport = getMappedCoverage(v,dictsam)
						v.addNormalSupport(normsupport)
					#	v.setUniqMapp(uniq)
					#	v.setOverallMapp(total)
						validmut.append(v)
	return validmut
def mutationSam(dictsam,genomes):
	try:
		mutations  = {}
		for k,v in dictsam.items():
			for align in v:
				muts = checkAlign(align,genomes)
				for e in muts:
		#			e.addNormalSupport(getMappedCoverage(e,dictsam))
					if e.__hash__() in mutations:
						mutations[e.__hash__()].addSupport(e.getSupport())
						#mutations[e.__hash__()].addNormalSupport(e.getNormalSupport())
				#		for s in e.getSupport():
				#			if s not in mutations[e.__hash__()].getSupport():
				#				mutations[e.__hash__()].addSupport(s)
					else:
						mutations[e.__hash__()] = e
		return mutations				
	except Exception as err:
		print(err)
		sys.exit(0)
def checkAlign(ali,genome):
	try:
		cigar = ali.getCigar()
		seq = ali.getSeq()
		qual = len(ali.getSeq())*'#'
		mutations = []
		pos = ali.getPos()
		chr = ali.getChr()
		id = ali.getID()
		if 'D' in cigar.upper() or 'I' in cigar.upper():
			match = align_candidate(read(id,seq,qual))
			#ali = alignment(id,chr,orient,pos,seq,cigar)
			match.addmatch([(ali,'')])
			cut = re.split('M|I|D',cigar)
			cut.pop()
 			off = 0
			roff = 0
			moff = 0
			for p in cut:
				off = off + len(p)
				if cigar[off] == 'M':
					roff = roff + int(p)
					moff = moff + int(p)
					off = off + 1
				elif cigar[off] == 'I':
					mu = mutation(1,chr,pos+roff,seq[roff:roff+int(p)],[match])
					mutations.append(mu)
					roff = roff + int(p)
					off = off + 1
				elif cigar[off] == 'D':
					#print(chr)
					#print(genome[chr],len(genome[chr]))
					#print(genome[chr][0][pos+moff-1],pos + moff-1+int(p))
					mu = mutation(2,chr,pos+moff,genome[chr][0][pos+moff-1:pos+moff-1+int(p)],[match])
					mutations.append(mu)
					off = off + 1
					moff = moff + int(p)
		return mutations
	except Exception as err:
		print(err)
		sys.exit(-1)
def codingRegionMutation(mutations,outdir,refgene,dictsam):
	try:
		mutfile = ''
		suppfile = ''
		validmut = []
		for m in mutations:
			if refgene.inCodingRegion(m.mutationPos()):
				if m.NumOfSupportReads() >= 2:
					#mut.write('%s\t%s\t%s\t%s\t%s\t%s\t'%(m.mutationType(),m.mutationChr().split('ex')[0],m.mutationPos(),m.mutationSeq(),m.NumOfSupportReads(),refgene.genename()))
					m.setChromosome(refgene.chr())
					m.setChrMutRange(refgene)
					#m.setChrPos(refgene.getChrPos(m.mutationStartPos()))
					m.setGene(refgene.genename())
					uniq,total = getMappedCoverage(m,dictsam)
					m.setUniqMapp(uniq)
					m.setOverallMapp(total)
					mutfile = mutfile + m.__str__()+'\n'
					#rds = ','.join(map(lambda x: x.split('\t')[0],m.getSupport()))
					for s in m.getSupport():
						suppfile = suppfile + s.__rec__() + '\n'
					validmut.append(m)
		if mutfile:
			assert(suppfile)
			mut = open(os.path.join(outdir,m.mutationChr()+'.mut'),'w')
			supp = open(os.path.join(outdir,m.mutationChr()+'.supp'),'w')
			mut.write(mutfile)
			supp.write(suppfile)
			mut.close()
			supp.close() 
		return validmut
	except Exception as err:
		print(err)
		sys.exit(0)
def extractCodingMutation(refgene,inputdir,outdir,dictsam):
	try:
		rna = refgene.rnaname()
		valid = []
		#if os.path.exists(os.path.join(inputdir,rna + '.mut')) and not os.path.exists(os.path.join(outdir,rna+'.mut')):
		mutations = parseMutation(os.path.join(inputdir,rna+'.mut'),os.path.join(inputdir,rna+'.supp'))
		valid =  codingRegionMutation(mutations,outdir,refgene,dictsam)
		return valid
	except Exception as err:
		print(err)
		sys.exit(0)
def extractCodingMutDrive(refgenefile,mutationdir,outdir,origsamfile):
	try:
		refgenes = parseRefGene(refgenefile)
		dictsam = loadSamByPos(origsamfile)
		overall = []
		for k,rnas in refgenes.items():
			for rna in rnas:
				if os.path.exists(os.path.join(mutationdir,rna.rnaname() + '.mut')):
					overall = overall + extractCodingMutation(rna,mutationdir,outdir,dictsam)
		#print(len(overall))	
		overall.sort(reverse=True)
		mf = open(os.path.join(outdir,'overall.mut'),'w')
		ms = open(os.path.join(outdir,'overall.supp'),'w')	
		for m in overall:
			mf.write(m.__str__()+'\n')
			for s in m.getSupport():
				ms.write(s.__rec__() + '\n')
		mf.close()
		ms.close()
	except Exception as err:
		print(err)
		sys.exit(0)
def combineExtractCodingRegion(refgenefile,mutationdir,prefix,outdir,cutoff):
	try:
		refgenes = parseRefGene(refgenefile)
		overall = []
		mutations = parseMutation(os.path.join(mutationdir,prefix+'.mut'),os.path.join(mutationdir,prefix+'.supp'))
		for m in mutations:
			if not m.mutationChr().split('ex')[0] in refgenes:
				raise('can not find mrna id %s in refgene file'%(m.rnaname()))
			for e in refgenes[m.mutationChr().split('ex')[0]]:
				if e.rnaname() == m.mutationChr():
					if e.inCodingRegion(m.mutationPos()) and m.NumOfSupportReads() >= cutoff:
						m.setChromosome(e.getChr())
						m.setChrMutRange(e)
						#m.setChrPos(e.getChrPos(m.mutationStartPos()))
						m.setGene(e.genename())
						overall.append(m)

						
		if mutations:
			muth = open(os.path.join(outdir,prefix + '.csv','w'))
			supph = open(os.path.join(outdir,prefix+'.csv','w'))
			mutations = list(set(overall))
			mutations.sort(reverse=True)
			for e in mutations:
				muth.write(e.__str__()+'\n')
				for s in e.getSupport():
					supph.write(s.__rec__()+'\n')
			muth.close()
			supph.close()
	except Exception as err:
		print(err)
		sys.exit(0)
						
				
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
