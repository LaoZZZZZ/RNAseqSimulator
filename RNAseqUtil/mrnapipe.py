#!/usr/bin/python

import os
import sys

from mapping import mapping
from mapping import pairmapping
from read_func import splitReadsPair
from read_func import splitreads
from mutation_detect import pairdrive
from mutation_detect import loadfqimpl
from parseRefGene import loadRefMrna,parseRefGene
from validateMutation import combineExtractCodingRegion as extractCoding
from removeDuplicate import removeDuplicate
from annotation import annotateCsv
from multiprocessing import Process
from extract_mapped import source_parsing
from scan_reads import singleSample as filterUnmapped
from mapping import mapping
from filter_mapped_trunc import filterMapped
from sort_sam import sort_bam
from filter_original import filterOrigin
from kmerCount import kmerCount
import subprocess 
from myemail import email2me
from refGeneStruct import refGene
from parseRefGene import parseRefGene
from samDrive import samDriveOnline
import random
# extract all mapped reads that are mapped to the specified resion
def extract_mapped_reads(inputdir,target,region,outputdir):
	try:
		if os.path.isdir(inputdir):
			targetdir = os.path.join(inputdir,target)
			for f in os.listdir(targetdir):
				if os.path.isfile(os.path.join(targetdir,f)) and f.endswith('bam'):
					input = os.path.join(targetdir,f)
					output = os.path.join(outputdir,target + '.mapped')
					command = 'samtools view ' + input + ' ' + region + ' > ' + output
					print(command)
					pipe = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
					pipe.communicate()
					if pipe.returncode:
						print('fail to extract mapped reads, something wrong')
						raise Exception('error')
	except Exception as error:
		print(error)
		sys.exit(0)
# transform the original bam file to sam file
def bam2sam(inputdir,target,outputdir):
	try:
		if os.path.isdir(inputdir):
			targetdir = os.path.join(inputdir,target)
			for f in os.listdir(targetdir):
				if os.path.isfile(os.path.join(targetdir,f)) and f.endswith('bam'):
					input = os.path.join(targetdir,f)
					output = os.path.join(outputdir,target +  '.sam')
					command = 'samtools view ' + input +  ' > ' + output
					print(command)
					pipe = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
					pipe.communicate()
					if pipe.returncode:
						print('fail to extract mapped reads, something wrong')
						raise Exception('error')
	except Exception as error:
		print(error)
		sys.exit(0)	


# pipeline for a single sample
def pipeline(inputdir,target,region,genome,outputdir):
	# extract already mapped reads from the bam file
	outputmapped = os.path.join(outputdir,'mapped')	
	#extract_mapped_reads(inputdir,target,region,outputmapped)
	# transform the bam fiel to samfile for later use
	outputsam = os.path.join(outputdir,'samfile')
	#bam2sam(inputdir,target,outputsam)
	# filter out those unmapped reads and split it into two pieces
	outputunmapped = os.path.join(outputdir,'unmapped')
	filterUnmapped(inputdir,outputsam,target,outputunmapped,20,1)
	# map the original unmapped reads to reference genome
	outputmapping = os.path.join(outputdir,'mappingResult')	
	mapping(outputunmapped,target,genome,outputmapping)	
	# filter out those mapped reads pieces
	outputfilter = os.path.join(outputdir,'region_mapped')
	filterMapped(outputmapping,target,outputfilter)
	# sort the sam file and generate corresponding sorted bam file
	sort_bam(outputfilter,target)
		
	# find the original complete reads from the sam file
	outputkmer = os.path.join(outputdir,'kmer')
	filterOrigin(outputunmapped,outputfilter,target,outputkmer)
	# generate all kmers from those original reads that mapped to this region
	kmerCount(outputkmer,target,23,outputkmer)
	subject = target + ' pipeline finished'
	body = 'Please check the attachment'
	email2me(subject,body,['pipeline_log'])	
def rnaPipeline(inputdir,target,genome,outputdir):
	# map the original unmapped reads to reference genome
	outputmapping = os.path.join(outputdir,'mappingResult')	
	mapping(outputunmapped,target,genome,outputmapping)	
	# filter out those mapped reads pieces
	outputfilter = os.path.join(outputdir,'region_mapped')
	filterMapped(outputmapping,target,outputfilter)
	# sort the sam file and generate corresponding sorted bam file
	sort_bam(outputfilter,target)
		
	# find the original complete reads from the sam file
	outputkmer = os.path.join(outputdir,'kmer')
	filterOrigin(outputunmapped,outputfilter,target,outputkmer)
	# generate all kmers from those original reads that mapped to this region
	kmerCount(outputkmer,target,23,outputkmer)
	subject = target + ' pipeline finished'
	body = 'Please check the attachment'
	email2me(subject,body,['pipeline_log'])
		
def append2ref(rna,seq,outputh):
	try:
		outputh.write('>'+rna+'\n')
		outputh.write(seq)
		outputh.write('\n')
	except Exception as err:
		print(err)
		import sys
		sys.exit(0)
def traverseMrna(dir):
	rnas = [] 
	for f in os.listdir(dir):
		rnas.append(f.split('.')[0])
	return rnas			
	

def findMapped(genome,refGene,genomefile,outmapping,inputfqdir,target):
	from parseRefGene import refGenerator
	from parseRefGene import ref2file
	from buildindex import buildindex
	emptyf = os.path.join(genome,'emptyRNAseq.txt')
	inputref = os.path.join(genome,'ref')
	inputrefh = open(os.path.join(inputref,'mrna.fa'),'w')
	emptyref = open(emptyf,'w')
	outputindex = os.path.join(genome,'index')
	total = 0
	for rna in refGenerator(genome,refGene,genomefile):
		if not len(rna['seq']):
			emptyref.write(rna['gene'].__str__() + '\n')
		else:
			total = total + 1
			# write the coding region to a fa file named with the gene name
			append2ref(rna['gene'].rnaname(),rna['seq'],inputrefh)
			ref2file(rna['gene'].rnaname(),rna['seq'],inputref)
			buildindex(inputref,rna['gene'].rnaname(),outputindex)	
	inputrefh.close()
	print(total)
	# build the bowtie2 index for the coding region
	if not os.path.isfile(os.path.join(outputindex,'mrna.1.bt2')) or not os.path.isfile(os.path.join(outputindex,'mrna.2.bt2')):
		buildindex(inputref,'mrna',outputindex)
	# mapping the originally unmapped reads to the index
	mrnaref = os.path.join(outputindex,'mrna')
	mapping(inputfqdir,target,mrnaref,outmapping,'mrna',inputfqdir)
	#os.rename(os.path.join(inputfqdir,target+'.unmapped'),os.path.join(inputfqdir,target+'.fq'))
# maping those finally splitted  unmapped reads into each mrna, keep those mapped parts
def partialMapping(inputindexdir,outputdir,inputfqdir,target):
	rnas = traverseMrna(inputindexdir)
	for e in rnas:
		genome = os.path.join(inputindexdir,e)
		mapping(inputfqdir,target+'_unmapped_split',genome,outputdir,e,'',False)
def validateMutation(refdir,originalfq,samdir,mutationdir,rna,cutoff):
	from mutation_detect import drive
	reffile = os.path.join(refdir,rna + '.fa')
	samfile = os.path.join(samdir,rna + '.sam')
	if not os.path.exists(os.path.join(mutationdir,rna+'.mut')): 
		drive(originalfq,samfile,reffile,0,float('inf'),mutationdir)			
def findMutation(inputdir,originalfq,splitmappeddir,refdir,mutationdir,cutoff):
	rnas = traverseMrna(os.path.join(inputdir,splitmappeddir))
	abrefdir = os.path.join(inputdir,refdir)
	abssamdir = os.path.join(inputdir,splitmappeddir)
	absmutationdir = os.path.join(inputdir,mutationdir)
	print(abrefdir,abssamdir,absmutationdir)
	assert(os.path.isdir(abssamdir) and os.path.isdir(absmutationdir) and os.path.isdir(absmutationdir))
	for rna in rnas:
		validateMutation(abrefdir,originalfq,abssamdir,absmutationdir,rna,cutoff)
def singleRNA(inputindexdir,mappingoutdir,inputfqdir,originalfq,splitfq,refdir,mutationdir,rna,cutoff):
	if not os.path.exists(os.path.join(mappingoutdir,rna+'.sam')):	
		from mapping import mapping
		genome = os.path.join(inputindexdir,rna)
		mapping(inputfqdir,splitfq,genome,mappingoutdir,rna,'',False)
	else:
		print(rna)
	import time
	time.sleep(1)	
	while not os.path.exists(os.path.join(mappingoutdir,rna+'.sam')):
		time.sleep(4)
	validateMutation(refdir,originalfq,mappingoutdir,mutationdir,rna,cutoff)
				
def batchRNA(workdir,inputindexdir,mappingoutdir,inputfqdir,originalfq,splitfq,refdir,mutationdir,cdmutationdir,rnas,cutoff):
	
	print(inputindexdir,mappingoutdir,inputfqdir,splitfq,refdir,mutationdir)
	total = 0
	refgenes = parseRefGene(os.path.join(workdir,'chromo/refGene.txt'))
	for e in rnas:
		total = total + 1
#		print(inputindexdir,mappingoutdir,inputfqdir,originalfq,splitfq,refdir,mutationdir,rnas)
		singleRNA(inputindexdir,mappingoutdir,inputfqdir,originalfq,splitfq,refdir,mutationdir,e,cutoff)
		
		assert(e.split('ex')[0] in refgenes)
		for mrna  in refgenes[e.split('ex')[0]]:
			if mrna.rnaname() == e:
				extractCodingMutation(mrna,mutationdir,cdmutationdir)
		if total %50 == 0:
			subject ='finished 50 mrna '
			body = 'Please check the attachment'
			email2me(subject,body)
def combineMapping(genome,refGene,genomefile,outmapping,inputfqdir,targets):
	from parseRefGene import refGenerator
	from parseRefGene import ref2file
	from buildindex import buildindex
	emptyf = os.path.join(genome,'emptyRNAseq.txt')
	inputref = os.path.join(genome,'ref')
	emptyref = open(emptyf,'w')
	outputindex = os.path.join(genome,'index')
	total = 0
	if not os.path.isfile(os.path.join(inputref,'mrna.fa')):
		inputrefh = open(os.path.join(inputref,'mrna.fa'),'w')
		for rna in refGenerator(genome,refGene,genomefile):
			if not len(rna['seq']):
				emptyref.write(rna['gene'].__str__() + '\n')
			else:
				total = total + 1
				# write the coding region to a fa file named with the gene name
				append2ref(rna['gene'].rnaname(),rna['seq'],inputrefh)
		inputrefh.close()
	# build the bowtie2 index for the coding region
	if not os.path.isfile(os.path.join(outputindex,'mrna.1.bt2')) or not os.path.isfile(os.path.join(outputindex,'mrna.2.bt2')):
		buildindex(inputref,'mrna',outputindex)
	# mapping the originally unmapped reads to the index
	mrnaref = os.path.join(outputindex,'mrna')
	for target in targets:
		mapping(inputfqdir,target,mrnaref,outmapping,target,inputfqdir)
	#os.rename(os.path.join(inputfqdir,target+'.unmapped'),os.path.join(inputfqdir,target+'.fq'))

def updatedSingleSampleProcessing(inputdir,target,genome,refgenes,refmrna,outputdir,sublen,split,nproc,cutoff):
	if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target + '_unique.csv')):
		inputfqdir = os.path.join(inputdir ,'unmapped')
		#bamfile = os.path.join(bamfiledir,target + '/' + target + '_accepted_hits.bam')
		if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target + '_small.csv')):
		#	if not os.path.isfile(os.path.join(inputdir,'mapped/' + target + '.sam')):
		#		mapping(inputfqdir,target,os.path.join(inputdir,'chromo/index/mrna'),os.path.join(inputdir,'mapped'),target,inputfqdir)
			samDriveOnline(os.path.join(inputdir,'mapped/'+target+'.sam'),refmrna,os.path.join(inputdir,'finalmutation'),target + '_small',cutoff,refgenes,nproc,split)	
		if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target + '_small_unique.csv')):
			#annotateCsv(os.path.join(inputdir,'finalmutation/' + target + '_small.csv'),refgenes,os.path.join(inputdir,'mapped/' + target + '.sam'),os.path.join(inputdir,'finalmutation/' + target + '_small_unique_annotated.csv'))	
			removeDuplicate(os.path.join(inputdir,'finalmutation/' + target + '_small.csv'),os.path.join(inputdir,'finalmutation/' + target + '_small.supp'),os.path.join(inputdir,'finalmutation/' + target + '_small.norm'))
	#	fqfile = os.path.join(inputfqdir, target+'_unmapped.fq')
	#	fq1 = target+'_split_1.fq'
	#	fq2 = target+'_split_2.fq'
		#output = os.path.join(inputfqdir,target+'_split.fq')
	#	if not os.path.isfile(os.path.join(inputfqdir,fq1)) or not os.path.isfile(os.path.join(inputfqdir,fq2)):
	#		splitReadsPair(fqfile,inputfqdir,target + '_split',sublen)
	#	if not os.path.isfile(os.path.join(outputdir,target+'_split.sam')):
	#		pairmapping(inputfqdir,fq1,fq2,os.path.join(genome,'index/mrna'),outputdir,target+'_split')
	#	splitsamfile = os.path.join(outputdir,target+'_split.sam')
		#print(inputdir,bamfiledir,target)
	#	dictfq = loadfqimpl(os.path.join(inputfqdir,target+'.fq'))
	#	if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target+'_large.csv')):
	#		pairdrive(dictfq,splitsamfile,origsamfile,refmrna,refgenes,os.path.join(inputdir,'finalmutation'),cutoff)
	#		os.system('cp ' + os.path.join(inputdir,'finalmutation/' + target+'.csv')  + '  ' + os.path.join(inputdir,'finalmutation/' + target+'_large.csv') )	
	#	if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target+'_large_unique_annotated.csv')):
	#		removeDuplicate(os.path.join(inputdir,'finalmutation/' + target + '_large.csv'))
	#		annotateCsv(os.path.join(inputdir,'finalmutation/' + target + '_large_unique.csv'),refgenes,bamfile,os.path.join(inputdir,'finalmutation/' + target + '_large_unique_annotated.csv'))	
	#	os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_small.csv') + ' >>  ' + os.path.join(inputdir,'finalmutation/' + target+'.csv')) 
	#	os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_large_unique.csv') + ' >  ' + os.path.join(inputdir,'finalmutation/' + target+'_unique.csv')) 
	#	os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_small_unique.csv') + ' >>  ' + os.path.join(inputdir,'finalmutation/' + target+'_unique.csv')) 
	#	os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_large_unique_annotated.csv') + ' >  ' + os.path.join(inputdir,'finalmutation/' + target+'_unique_annotated.csv')) 
	#	os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_large_unique_annotated.csv') + ' >>  ' + os.path.join(inputdir,'finalmutation/' + target+'_unique_annotated.csv')) 
		#removeDuplicate(os.path.join(inputdir,'finalmutation/' + target + '.csv'))
		#annotateCsv(os.path.join(inputdir,'finalmutation/' + target + '_unique.csv'),refgenes,bamfile,os.path.join(inputdir,'finalmutation/' + target + '_unique_annotated.csv'))	
	subject = target + ' mrnapipeline finished'
	body = 'Please check the attachment'
	email2me(subject,body)	
def singleSampleProcessing(inputdir,target,genome,refgenes,refmrna,outputdir,sublen,bamfiledir,cutoff):
	if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target + '_unique.csv')):
		inputfqdir = os.path.join(inputdir ,'unmapped')
		origf = open(os.path.join(inputdir,target + '.sam'),'w')
		origf.write('@' + target + '\n')
		origf.close()
		origsamfile = os.path.join(inputdir,target + '.sam')
		bamfile = os.path.join(bamfiledir,target + '/' + target + '_accepted_hits.bam')
		if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target + '_small.csv')):
			if not os.path.isfile(os.path.join(inputdir,'mapped/' + target + '.sam')):
				mapping(inputfqdir,target,os.path.join(inputdir,'chromo/index/mrna'),os.path.join(inputdir,'mapped'),target,inputfqdir)
			samDriveOnline(os.path.join(inputdir,'mapped/'+target+'.sam'),refmrna,origsamfile,os.path.join(inputdir,'finalmutation'),target + '_small',cutoff,refgenes)	
		if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target + '_small_unique_annotated.csv')):
			removeDuplicate(os.path.join(inputdir,'finalmutation/' + target + '_small.csv'))
			annotateCsv(os.path.join(inputdir,'finalmutation/' + target + '_small_unique.csv'),refgenes,bamfile,os.path.join(inputdir,'finalmutation/' + target + '_small_unique_annotated.csv'))	
		fqfile = os.path.join(inputfqdir, target+'_unmapped.fq')
		fq1 = target+'_split_1.fq'
		fq2 = target+'_split_2.fq'
		#output = os.path.join(inputfqdir,target+'_split.fq')
		if not os.path.isfile(os.path.join(inputfqdir,fq1)) or not os.path.isfile(os.path.join(inputfqdir,fq2)):
			splitReadsPair(fqfile,inputfqdir,target + '_split',sublen)
		if not os.path.isfile(os.path.join(outputdir,target+'_split.sam')):
			pairmapping(inputfqdir,fq1,fq2,os.path.join(genome,'index/mrna'),outputdir,target+'_split')
		splitsamfile = os.path.join(outputdir,target+'_split.sam')
		#print(inputdir,bamfiledir,target)
		dictfq = loadfqimpl(os.path.join(inputfqdir,target+'.fq'))
		if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target+'_large.csv')):
			pairdrive(dictfq,splitsamfile,origsamfile,refmrna,refgenes,os.path.join(inputdir,'finalmutation'),cutoff)
			os.system('cp ' + os.path.join(inputdir,'finalmutation/' + target+'.csv')  + '  ' + os.path.join(inputdir,'finalmutation/' + target+'_large.csv') )	
		if not os.path.isfile(os.path.join(inputdir,'finalmutation/' + target+'_large_unique_annotated.csv')):
			removeDuplicate(os.path.join(inputdir,'finalmutation/' + target + '_large.csv'))
			annotateCsv(os.path.join(inputdir,'finalmutation/' + target + '_large_unique.csv'),refgenes,bamfile,os.path.join(inputdir,'finalmutation/' + target + '_large_unique_annotated.csv'))	
		os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_small.csv') + ' >>  ' + os.path.join(inputdir,'finalmutation/' + target+'.csv')) 
		os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_large_unique.csv') + ' >  ' + os.path.join(inputdir,'finalmutation/' + target+'_unique.csv')) 
		os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_small_unique.csv') + ' >>  ' + os.path.join(inputdir,'finalmutation/' + target+'_unique.csv')) 
		os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_large_unique_annotated.csv') + ' >  ' + os.path.join(inputdir,'finalmutation/' + target+'_unique_annotated.csv')) 
		os.system('cat ' + os.path.join(inputdir,'finalmutation/' + target + '_large_unique_annotated.csv') + ' >>  ' + os.path.join(inputdir,'finalmutation/' + target+'_unique_annotated.csv')) 
		#removeDuplicate(os.path.join(inputdir,'finalmutation/' + target + '.csv'))
		#annotateCsv(os.path.join(inputdir,'finalmutation/' + target + '_unique.csv'),refgenes,bamfile,os.path.join(inputdir,'finalmutation/' + target + '_unique_annotated.csv'))	
	subject = target + ' mrnapipeline finished'
	body = 'Please check the attachment'
	email2me(subject,body)	
# first map the whole reads into the rna index
# second map the split reads to the rna index
# identify possible mutation using the split mapping
def rnapipeCombine(inputdir,targets,genome,outputdir,sublen,bamfiledir,nsplit,nproc,cutoff = 2):
	from refGeneStruct import refGene
	import random
	# the directory that hold those mapped sam file for each mrna sequence
	outmapping = os.path.join(inputdir,'mapped')
	# output those unmapped reads into this directory
	refmrna = loadRefMrna(os.path.join(genome,'ref/mrna.fa'))
	refgenes = parseRefGene(os.path.join(genome,'refGene.txt'))
	size = int(len(refgenes)/nsplit)
	left = len(refgenes)%nsplit
	split = set()
	for k, v in refgenes.items():
		for r in v:
			split.add(r.rnaname())
	from samDrive import splitMrna
	#splits = splitMrna(split,nsplit)
	# filter out those mapped reads
	#combineMapping(genome,'refGene.txt','genome.fa',outmapping,inputfqdir,targets)
	# split the those unmapped reads into three pieces, front, end and middle, only the middle part will be discarded
	proc = []
	from multiprocessing import Process
	for target in targets:
		print(target)
		p = Process(target=updatedSingleSampleProcessing,args=((inputdir,target,genome,refgenes,refmrna,outputdir,sublen,split,nproc,cutoff)))	
		p.start()
		proc.append(p)
		if len(proc) == 4:
			for e in proc:
				e.join()
			proc = []
	if proc:
		for e in proc:
			e.join()
	subject = ' mrnapipeline finished'
	body = 'Please check the attachment'
	email2me(subject,body)	
def wholeRNA(inputdir,targets,genome,outputdir,sublen,nproc,cutoff = 2):
	from mapping import mapping
	import os
	import sys
	# the directory that hold those mapped sam file for each mrna sequence
	outmapping = os.path.join(inputdir,'mapped')
	# output those unmapped reads into this directory
	inputfqdir = os.path.join(inputdir ,'unmapped')
	# filter out those mapped reads
	findMapped(genome,'refGene.txt','genome.fa',outmapping,inputfqdir,targets[0])
	# split the those unmapped reads into three pieces, front, end and middle, only the middle part will be discarded
	from read_func import splitreads
	fqfile = os.path.join(inputfqdir,targets[0]+'_unmapped.fq')
	output = os.path.join(inputfqdir,targets[0]+'_unmapped_split.fq')
	splitreads(fqfile,output,sublen)
	rnas = traverseMrna(os.path.join(inputdir,'chromo/ref'))
	rnas.remove('mrna')
	exitrnas = traverseMrna(os.path.join(inputdir,'splitMapped'))
	pending = [val for val in rnas if not val in exitrnas]	
	print(len(exitrnas),len(rnas),len(pending))	
	size = len(pending)/nproc
	left = len(pending)%nproc
	split = []
	for i in range(nproc):
		split.append(pending[i*size:(i+1)*size])
	split[-1] = exitrnas + split[-1]
	if left:
		split[-1] = split[-1] + pending[nproc*size:]
	from mutation_detect import loadfqimpl
	import os
	originalfq = loadfqimpl(fqfile)
	mappingdir = os.path.join(inputdir,'splitMapped')
	splitfq = os.path.basename(output).split('.')[0]
	refdir = os.path.join(inputdir,'chromo/ref')
	mutationdir = os.path.join(inputdir,'mutation')
	inputindexdir = os.path.join(inputdir,'chromo/index')
	cdmutationdir = os.path.join(inputdir,'cdmutation')
	proc = []
	from multiprocessing import Process
	for i in range(len(split)):
		p = Process(target=batchRNA,args=((inputdir,inputindexdir,mappingdir,inputfqdir,originalfq,splitfq,refdir,mutationdir,cdmutationdir,split[i],cutoff)))	
		p.start()
		proc.append(p)
	# find those partially mapped reads for each mrna sequence
	#partialMapping(os.path.join(inputdir,'chromo/index'),os.path.join(inputdir,'splitMapped'),inputfqdir,targets[0])
	# validate those alignment and find possible mutation for each mrnasequence	
	#findMutation(inputdir,os.path.join(inputfqdir,targets[0]+'_unmapped.fq',splitMapped','chromo/ref','mutation')
	for e in proc:
		e.join()	
	subject = targets[0] + ' mrnapipeline finished'
	body = 'Please check the attachment'
	email2me(subject,body)	
if __name__ == '__main__':
	task = source_parsing(sys.stdin)
	outputdir = '/home/luzhao/THR104_analysis'
	print(task['dir'],task['targets'],task['genome'])
	sublen = 20
	inputdir = task['dir']
	targets = task['targets']
	genome = task['genome']
	inputfqdir = os.path.join(inputdir,'unmapped')
	#fqfile = os.path.join(inputfqdir,target+'.fq')
	#output = os.path.join(inputfqdir,target+'_unmapped__split.fq')
	#splitreads(fqfile,output,sublen)
	#findMutation(inputdir,'splitMapped','chromo/ref','mutation')
	#partialMapping(os.path.join(inputdir,'chromo/index'),os.path.join(inputdir,'splitMapped'),inputfqdir,target)
	cutoff = 2
	bamdir = '/home/stc/Documents/GenomicData/RNAseq/WBahou'
	nproc = 1
	nsplit = 1
	print(targets)
	rnapipeCombine(inputdir,targets,genome,os.path.join(outputdir,'partialMapping'),sublen,bamdir,nsplit,nproc,cutoff)
#	wholeRNA(task['dir'],task['targets'],task['genome'],outputdir,sublen,bamdir,cutoff)
	#for target in task['targets']:
	#	pipeline(task['dir'],target,region,task['genome'],outputdir)
	#	Process(target=pipeline,args=(task['dir'],target,region,task['genome'],outputdir)).start()
	#print('Main thread exit' + '\n')

