#!/usr/bin/python

from refGeneStruct import refGene
from parseRefGene import parseRefGene
from parseRefGene import parseMutation
from mutation import mutation
from align_candidate import align_candidate
from alignment import alignment
from strMatch import *
from reads import read
import os
import sys
import re
from parseRefGene import loadRefMrna
def samDrive(samfile,genomefile,outdir,prefix,cutoff):
	mutd = mutationSam(samfile,genomefile)
	mutations = filterMutations(mutd,cutoff)
	outputMutations(outdir,prefix,mutations)	
def outputMutations(outdir,prefix,mutations):
	mutfile = os.path.join(outdir,prefix+'.mut')
	suppfile = os.path.join(outdir,prefix+'.supp')	
	print(mutfile,suppfile)
	if mutations:
		muth = open(mutfile,'w')
		supph = open(suppfile,'w')
		for e in mutations:
			muth.write(e.__str__()+'\n')
			for s in e.getSupport():
				supph.write(s.__rec__()+'\n')
		muth.close()
		supph.close()
	
def filterMutations(mutations,cutoff):
	validmut = []
	if mutations:
		for k,v in mutations.items():
			if len(v.getSupport()) >= cutoff:
				validmut.append(v)
	return validmut
def mutationSam(file,genomefile):
	try:
		mutations = {}
		fhandler = open(file,'r')
		genomes = loadRefMrna(genomefile)	
		for line in fhandler:
			if line[0] != '@':
				align = line.rstrip('\n').split('\t')
				muts = checkAlign(align,genomes)
				for e in muts:
					if e.__hash__() in mutations:
						for s in e.getSupport():
							if s not in mutations[e.__hash__()].getSupport():
								mutations[e.__hash__()].addSupport(e.getSupport())
					else:
						mutations[e.__hash__()] = e
		return mutations				
	except Exception as err:
		print(err)
		sys.exit(0)
def checkAlign(align,genome):
	try:
		pos = int(align[3])
		id = align[0]
		cigar = align[5]
		seq = align[9]
		qual  = align[10]
		chr = align[2]
		orient = align[1]
		mutations = []
		#print(id,orient,chr,pos,cigar,seq,qual)
		if cigar.upper().count('D') or cigar.upper().count('I'):
			match = align_candidate(read(id,seq,qual))
			ali = alignment(id,chr,orient,pos,seq)
			match.addmatch([(ali,'')])
			cut = re.split('M|I|D',cigar)[0:-2]
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
def codingRegionMutation(mutations,outdir,refgene):
	try:
		for m in mutations:
			if refgene.inCodingRegion(m.mutationPos()):
				distinctrd = set() 
				for s in m.getSupport():
					rec = s.split()
					id = int(rec[3])
					distinctrd.add(id)
				if len(distinctrd) > 2:
					mut = open(os.path.join(outdir,m.mutationChr()+'.mut'),'w')
					mut.write('%s\t%s\t%s\t%s\t%s\t%s\t'%(m.mutationType(),m.mutationChr().split('ex')[0],m.mutationPos(),m.mutationSeq(),m.NumOfSupportReads(),refgene.genename()))
					rds = ','.join(map(lambda x: x.split('\t')[0],m.getSupport()))
					mut.write(rds+'\n')
					mut.close()
					supp = open(os.path.join(outdir,m.mutationChr()+'.supp'),'w')
					for s in m.getSupport():
						supp.write(s)
					supp.close() 

	except Exception as err:
		print(err)
		sys.exit(0)
def extractCodingMutation(refgene,inputdir,outdir):
	try:
		rna = refgene.rnaname()
		print(rna)
		if os.path.exists(os.path.join(inputdir,rna + '.mut')) and not os.path.exists(os.path.join(outdir,rna+'.mut')):
			mutations = parseMutation(os.path.join(inputdir,rna+'.mut'),os.path.join(inputdir,rna+'.supp'))
			codingRegionMutation(mutations,outdir,refgene)
	except Exception as err:
		print(err)
		sys.exit(0)
def extractCodingMutDrive(refgenefile,mutationdir,outdir):
	try:
		refgenes = parseRefGene(refgenefile)
		for k,rnas in refgenes.items():
			for rna in rnas:
				if os.path.exists(os.path.join(mutationdir,rna.rnaname() + '.mut')):
					extractCodingMutation(rna,mutationdir,outdir)
	except Exception as err:
		print(err)
		sys.exit(0)
		
if __name__ == '__main__':
	dir = '/home/luzhao/THR104_analysis'
	refgene = os.path.join(dir,'chromo/refGene.txt')
	inputdir = os.path.join(dir,'mutation')
	outdir = os.path.join(dir,'finalmutation')
	extractCodingMutDrive(refgene,inputdir,outdir)
	#ref = os.path.join(dir,'chromo/ref/mrna.fa')
	#inputdir = os.path.join(dir,'mapped')
	#samfile = os.path.join(dir,'mapped/mrna.sam')
	#samDrive(samfile,ref,inputdir,'mrna',2)
	#mutationSam(samfile,ref)
