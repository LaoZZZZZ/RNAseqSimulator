#!/usr/bin/python

import os
import sys

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
	bam2sam(inputdir,target,outputsam)
	# filter out those unmapped reads and split it into two pieces
	outputunmapped = os.path.join(outputdir,'unmapped')
	#filterUnmapped(inputdir,outputsam,target,outputunmapped,20,1)
	# map the original unmapped reads to reference genome
	outputmapping = os.path.join(outputdir,'mappingResult')	
	#mapping(outputunmapped,target,genome,outputmapping)	
	# filter out those mapped reads pieces
	outputfilter = os.path.join(outputdir,'region_mapped')
	#filterMapped(outputmapping,target,outputfilter)
	# sort the sam file and generate corresponding sorted bam file
	#sort_bam(outputfilter,target)
		
	# find the original complete reads from the sam file
	#outputkmer = os.path.join(outputdir,'kmer')
	#filterOrigin(outputunmapped,outputfilter,target,outputkmer)
	# generate all kmers from those original reads that mapped to this region
	#kmerCount(outputkmer,target,23,outputkmer)
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
def wholeRNA(inputdir,targets,genome,outputdir):
	from parseRefGene import refGenerator
	from buildindex import buildindex
	genome = os.path.join(inputdir,chromo)
	for rna in refGenerator(genome,'refGene.txt','refMrna.fa'):
		buildindex(genome,rna)
		rnaRef = genome + rna
		for target in targets:
			rnaPipeline(inputdir,target,rnaRef,outputdir) 
if __name__ == '__main__':
	chr = 'chr19'
	range = '13048983-13056060'
	region = chr + ':' + range
	task = source_parsing(sys.stdin)
	outputdir = '/home/luzhao/THR104_analysis'
	for target in task['targets']:
		pipeline(task['dir'],target,region,task['genome'],outputdir)
	#	Process(target=pipeline,args=(task['dir'],target,region,task['genome'],outputdir)).start()
	print('Main thread exit' + '\n')
	
