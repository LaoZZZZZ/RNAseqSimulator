#!/usr/bin/python

from parseMutSupp import parseMutSupp
def findDuplicate(file):
	import os
	import sys
	dictrna = {}
	total = 0
	try:
		if os.path.isfile(file):
			fhandler = open(file,'r')
			for line in fhandler:
				#print(line)
				total = total + 1
				rec = line.rstrip('\n').split()
				if rec[1] in dictrna:
					print(rec[1])
				else:
					dictrna[rec[1]] = 1
		else:
			raise Exception('invalid file name ' + file)
		return dictrna
	except Exception as err:
		print(err)
		sys.exit(-1)
def parseRefGeneByGene(file):
	import os
	import sys
	from refGeneStruct import refGene
	dictrna = {}
	total = 0
	try:
		if os.path.isfile(file):
			fhandler = open(file,'r')
			for line in fhandler:
				rec = line.rstrip('\n').split()
				total = total + 1
				if rec[12] in dictrna:
					dictrna[rec[12]].append(refGene(rec[1],rec[2],rec[4],rec[5],rec[6],rec[7],rec[3],rec[8],rec[9].rstrip(',').split(','),rec[10].rstrip(',').split(','),rec[12]))	
				else:
					dictrna[rec[12]] = [refGene(rec[1],rec[2],rec[4],rec[5],rec[6],rec[7],rec[3],rec[8],rec[9].rstrip(',').split(','),rec[10].rstrip(',').split(','),rec[12])]	
		else:
			raise Exception('invalid file name ' + file)
		return dictrna
	except Exception as err:
		print(err)
		sys.exit(-1)
def parseRefGene(file):
	import os
	import sys
	from refGeneStruct import refGene
	dictrna = {}
	total = 0
	try:
		if os.path.isfile(file):
			fhandler = open(file,'r')
			for line in fhandler:
				rec = line.rstrip('\n').split()
					
				total = total + 1
				if rec[1] in dictrna:
#					raise Exception('the gene id(%s) is not unique'%(rec[1]))
					dictrna[rec[1]].append(refGene(rec[1],rec[2],rec[4],rec[5],rec[6],rec[7],rec[3],rec[8],rec[9].rstrip(',').split(','),rec[10].rstrip(',').split(','),rec[12]))	
				else:
					dictrna[rec[1]] = [refGene(rec[1],rec[2],rec[4],rec[5],rec[6],rec[7],rec[3],rec[8],rec[9].rstrip(',').split(','),rec[10].rstrip(',').split(','),rec[12])]	
		else:
			raise Exception('invalid file name ' + file)
		return dictrna
	except Exception as err:
		print(err)
		sys.exit(-1)
def parseRefMrna(rna,file):
	try:
		find = False
		ref = ''
		for line in file:
			if line[0] == '>' and not find:
				rec = line.rstrip('\n').split()
				if rec[0][1:] == rna:
					find = True
			elif line[0]== '>' and find and not ref:
				raise Exception('invalid fa format')
			elif find and line[0] != '>':
				ref = ref + line.rstrip('\n')
			elif ref and find and line[0] == '>':
				break
			else:
				pass
		return ref
	except Exception as err:
		print(err)
		sys.exit(1)
def loadRefMrna(file):
	try:
		mrna = {}
		handler = open(file,'r')
		next = False
		id = ''
		ref = ''
		for line in handler:
			if line[0] == '>' and not next:
				rec = line.rstrip('\n').split()
				id = rec[0][1:]
				next = True
			elif line[0]== '>' and next and not ref:
				raise Exception('invalid fa format')
			elif next and line[0] != '>':
				ref = ref + line.rstrip('\n')
			elif ref and next and line[0] == '>':
				if id in mrna:
					mrna[id].append(ref)
				else:
					mrna[id]=[ref]
				id = line.rstrip('\n').split()[0][1:]
				next = True
				ref = ''	
			else:
				pass
		if ref and next:
			if id in mrna:
				mrna[id].append(ref)
			else:
				mrna[id] = [ref]
		return mrna
	except Exception as err:
		print(err)
		sys.exit(0)
def ref2file(rna,seq,outdir):
	import os
	outputfile = os.path.join(outdir,rna+'.fa')	
	out = open(outputfile,'a+')
	out.write('>'+rna+'\n')
	out.write(seq)
	out.close()	

def refGenerator(inputdir,refgene,refMrna):
	import os
	try:
		indexgene = os.path.join(inputdir,refgene)
		rna = parseRefGene(indexgene)
		mrna = os.path.join(inputdir,refMrna)
		mrnaf = loadRefMrna(mrna)
		print('from refGenerator %s'%(len(rna)))
		for k in rna.keys():
			for e in rna[k]:
				if not e.chr() in mrnaf:
					yield {'gene':e,'seq':''}	
					#raise Exception('the chromosome(%s) of %s is not in the mrna dictionary '%(rna[k].chr(),k))
				else:
					#stretch = e.getCodingRegion(mrnaf[e.chr()][0])
					# get the rna seq reference
					stretch = e.getRNASeq(mrnaf[e.chr()][0])
					yield {'gene':e,'seq':stretch}
	except Exception as err:
		print(err)
		import sys
		sys.exit(0)
	
def duplicates(inputdir,refgene,refMrna):
	indexgene = os.path.join(inputdir,refgene)
	rna = parseRefGene(indexgene)
	mrna = os.path.join(inputdir,refMrna)
	mrna = loadRefMrna(mrna)
	dup = 'duplicate.fa'
	out = open(os.path.join(inputdir,dup),'a+')
	for k in rna.keys():
		stretch = mrna[k]
		if len(rna[k]) > 1 and len(stretch) == 1:
			if stretch:
				out.write('>'+k+'\n')
				out.write(stretch[0]+'\n')
			print(k)
	out.close()			
	
if  __name__ == '__main__':
	direct = '/home/luzhao/THR104_analysis/chromo'
	file = 'refGene.txt'
	import os
	import sys
	from buildindex import buildindex
	#genome = loadRefMrna(os.path.join(direct,'chr2.fa'))
	print('running')
	#print(len(genome))
	#for k,v in genome.items():
	#	print(k,len(v))
	rna = parseRefGeneByGene(os.path.join(direct,'refGene.txt'))
	print(len(rna))
	total = 0
	for k,v in rna.items():
		total += len(v)
	print(total)	

	#for i in range(50):
	#	print(sr.getChrPos(i))	
	#mrna = loadRefMrna(os.path.join(direct,'refMrna.fa'))	
	#original = '/home/luzhao/THR104_analysis'
	#for k in rna.keys():
	#	if rna[k].chr() == 'chr2' and rna[k].orient() == '+':
	#		print(k)
	#		os = open(os.path.join(original,k+'.s'),'w')
	#		os.write('>'+k+'\n')
	#		os.write(rna[k].getCodingRegion(genome[rna[k].chr()][0]).upper().rstrip('\n'))
	#		os.close()
	#		break			
	#duplicates(direct,file,'refMrna.fa')

