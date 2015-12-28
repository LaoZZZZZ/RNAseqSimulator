#!/usr/bin/python
import os
import sys
from parseMutSupp import parseMutation
from parseRefGene import parseRefGene,parseRefGeneByGene
from mutation import mutation
from mutation_detect import loadSamByPos
from mutation_detect import getMappedCoverage,loadPartialSamByPos
from extract_mapped import extractRegionMapping
from removeDuplicate import removeDuplicate
def countMapped(bamfile,region,outfile):
	extractRegionMapping(bamfile,region,outfile)
	dictsam = loadSamByPos(outfile)
	total = 0
	for k,v in dictsam.items():
		total = total + v
	return (len(dictsam),total)	
def annotate(mutfile,refgenes,orig_bamfile,outfile):
	muthandler = open(mutfile,'r')
	#dictsam = loadSamByPos(orig_samfile)
	newrec = ''
	out = open(outfile,'w')
	basename = os.path.splitext(os.path.basename(mutfile))[0]
	for line in muthandler:
		rec = line.rstrip('\n').split(',')
		newrec = rec[0] + ',' + rec[1] + ','+rec[2]+','+rec[3]
		rna = rec[1]
		rnapos = int(rec[2])
		muttype = int(rec[0])
		mutseq = rec[6]
		uniqsupp = int(rec[7])
		totalsupp = int(rec[8])
		mut = mutation(muttype,rna,rnapos,mutseq,[])
		if rna.split('ex')[0] in refgenes:
			for e in refgenes[rna.split('ex')[0]]:
				if e.rnaname() == rna:
					mut.setChrMutRange(e)
					start,end = mut.getChrMutRange()
					newrec = newrec + ',' + str(start) + ',' + str(end) + ',' + e.genename() + ',' + str(len(mutseq)) + ',' + mutseq + ','
					mut.setChromosome(e.chr())
					mut.setGene(e.genename())
				#	if orig_bamfile:
					#uniq,total = getMappedCoverage(mut,dictsam)
					uniq,total = countMapped(orig_bamfile,e.chr()+':' + str(mut.getChrMutRange()[0])+'-'+str(mut.getChrMutRange()[1]),basename)
					newrec = newrec + str(uniqsupp) + ',' + str(totalsupp) + ',' + str(uniq) + ',' + str(total) + ','
					newrec = newrec + str(float(uniqsupp)/(uniq + uniqsupp)) + ',' + str(float(totalsupp)/(totalsupp + total)) + ',' +  rec[-1] + '\n'
					out.write(newrec)
	os.system('rm ' + basename)
	out.close()

def annotateCsv(mutfile,refgenes,orig_bamfile,outfile):
	muthandler = open(mutfile,'r')
	if orig_bamfile.split('.')[1] == 'sam':
		dictsam = loadSamByPos(orig_bamfile)
	newrec = ''
	out = open(outfile,'w')
	basename = os.path.splitext(os.path.basename(mutfile))[0]
	for line in muthandler:
		rec = line.rstrip('\n').split(',')
		if not rec[0].isdigit():
			out.write(line)
			continue
		newrec = rec[0] + ',' + rec[1] + ','+rec[2]+','+rec[3]
		rna = rec[1]
		rnapos = int(rec[2])
		muttype = int(rec[0])
		mutseq = rec[8]
		uniqsupp = int(rec[9])
		totalsupp = int(rec[10])
		mut = mutation(muttype,rna,rnapos,mutseq,[])
		if rna.split('ex')[0] in refgenes:
			for e in refgenes[rna.split('ex')[0]]:
				if e.rnaname() == rna:
					mut.setChrMutRange(e)
					start,end = mut.getChrMutRange()
					newrec = newrec + ',' + str(start) + ',' + str(end) + ',' + e.genename() + ',' + str(len(mutseq)) + ',' + mutseq + ','
					mut.setChromosome(e.chr())
					mut.setGene(e.genename())
					if dictsam:
						uniq,total = getMappedCoverage(mut,dictsam)
					else:
						uniq,total = countMapped(orig_bamfile,e.chr()+':' + str(mut.getChrMutRange()[0])+'-'+str(mut.getChrMutRange()[1]),basename)
					newrec = newrec + str(uniqsupp) + ',' + str(totalsupp) + ',' + str(uniq) + ',' + str(total) + ','
					newrec = newrec + str(float(uniqsupp)/(uniq + uniqsupp)) + ',' + str(float(totalsupp)/(totalsupp + total)) + ',' +  rec[-1] + '\n'
					if len(newrec.split(',')) != 16:
						print(newrec,line)
						sys.exit(0)
					out.write(newrec)
					break
	os.system('rm '+ basename)
	out.close()
def findNormSupport(inputdir,prefix,samfile,refgenes):
	from mutation import mutation
	mutations = parseMutation(os.path.join(inputdir,prefix+'.csv'),os.path.join(inputdir,prefix+'.supp'),os.path.join(inputdir,prefix+'.norm'))
	mrnas = []
	for m in mutations:
		gene = m.getGene()
		tmp =  refgenes[gene]
		for r in tmp:
			mrnas.append(r.rnaname())
			print(r.rnaname())
	dictsam = loadPartialSamByPos(samfile,mrnas)
	normfile = open(os.path.join(inputdir,prefix+'.norm'),'w')
	for m in mutations:
		for mrna in refgenes[m.getGene()]:
			rna = mrna.rnaname()
			m.setRna(rna)
			m.setMutationPos(min(mrna.getMrnaPos(m.getChrMutRange()[0]),mrna.getMrnaPos(m.getChrMutRange()[1])))
			normsupports = getMappedCoverage(m,dictsam)
			for cand in normsupports:
				normfile.write(cand.__rec__() + '\n')
	normfile.close()
def postProcess(target):
	dir = '/home/luzhao/THR104_analysis'
	removeDuplicate(os.path.join(dir,'finalmutation/' + target + '.csv'))
	mutfile = os.path.join(dir,'finalmutation/' + target + '_unique.csv')
	refgenes = parseRefGene(os.path.join(dir,'chromo/refGene.txt'))
	origsamfile = '/home/stc/Documents/GenomicData/RNAseq/WBahou/' + target +  '/' + target + '_accepted_hits.bam'
	outfile = os.path.join(dir,'finalmutation/' + target + '_unique_anno.csv')
	annotateCsv(mutfile,refgenes,origsamfile,outfile)
if __name__ == '__main__':
	postProcess('THR104')
	dir = '/home/stc/Documents/GenomicData/RNAseq/WBahou/THR104'
	direct = '/home/luzhao/THR104_analysis/finalmutation'
	samfile = '/home/luzhao/THR104_analysis/mapped/THR102.sam'
	refgenes = parseRefGeneByGene('/home/luzhao/THR104_analysis/chromo/refGene.txt')
	findNormSupport(direct,'THR102_small_unique',samfile,refgenes)	
	#bamfile = os.path.join(dir,'THR104_accepted_hits.bam')
	#region = 'chr1:120612005-120612006'
	#outfile = '/home/luzhao/THR104_analysis/notch.sam'
	#countMapped(bamfile,region,outfile)
