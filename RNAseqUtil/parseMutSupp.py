#!/usr/bin/python

from alignment import alignment
from align_candidate import align_candidate
from reads import read
from mutation_detect import trim
from sampleInfo import sampleInfo
import os
import sys
from mutation import mutation
def parseMutation(mutfile,supportfile,normfile):
	try:
		mut = []
		print(mutfile,supportfile,normfile)
		if not os.path.isfile(mutfile) or not os.path.isfile(supportfile):
			raise Exception(' can not find files(%s,%s)'%(mutfile,supportfile))
		muthandler = open(mutfile,'r')
		mutations = []
		for line in muthandler:
			if line[0] == '@':
				continue
			rec = line.split(',')
			mutations.append(rec)
		mutsupports = parseMutSupp(supportfile)
		normals = parseMutSupp(normfile)
		for m in mutations:
			mutrds = m[15].split()
			normrds = m[16].split()
		#	print(m[16],normrds)
			mutmatches = []
			normmatches = []
			for r in mutrds:
				if not r in mutsupports:
					print(r)
					raise Exception('incompatible mutation file and mutation support file')
				mutmatches = mutmatches + mutsupports[r]
			for r in normrds:
				if not r in normals:
					print(r)
					raise Exception('incompatible  mutation  file and normal support file')
				normmatches = normmatches + normals[r]
					
			tmp = mutation(m[0],m[1],m[2],m[8],mutmatches)
			tmp.addNormalSupport(normmatches)
			tmp.setChrStart(int(m[4]))
			tmp.setChrEnd(int(m[5]))
			tmp.setChromosome(m[3])
			tmp.setGene(m[6])
			#tmp.setUniqMapp(int(m[11]))
			#tmp.setOverallMapp(int(m[12]))	
			mut.append(tmp)
#			for e in mut:
#				for a in e.getSupport():
#					print(a)
		return mut
	except Exception as err:
		print(err)
		sys.exit(0)
#parse the supp file, return a dictionary
def parseMutSupp(file):
	try:
		if not os.path.exists(file):
			raise Exception('can not find file(%s)' %(file))
		handler = open(file,'r')
		matches = {}
		cur = 0
		for line in handler:
			cur = cur + 1
			if line=='None\n' or line == 'na\n' or line == '\n':
				continue
			rec = line.rstrip('\n').split()
			if len(rec) != 5 and len(rec) != 6:
				print(len(rec),rec,cur)
				rec = line.split('HWI')
				print(rec)
				raise Exception(' %s seems not a mutation support file'%(file))
			if len(rec) == 5:
				al = alignment(rec[0],rec[1],rec[2],rec[3],rec[4])
			elif len(rec) == 6:
				al = alignment(rec[0],rec[1],rec[2],rec[3],rec[4],rec[5])
			if al.getID() in matches:
				matches[al.getID()].append(al)
			else:
				matches[al.getID()] = [al]
		res = {} 
		for k,v in matches.items():
			sv = sorted(v)
			seq = ''
			for al in sv:
				seq = seq + al.getSeq()	
			rd = read(v[0].getID(),seq,len(seq)*'a')
			m = align_candidate(rd)
			if len(sv) == 2:
				m.addmatch([tuple(sv)])
				res[v[0].getID()] = [m]
			elif len(sv) == 1:
				m.addmatch([tuple([sv[0], ''])])
				res[v[0].getID()] = [m]
			elif len(sv) > 2:
				als = []
				for e in sv:
					als.append(tuple([e,'']))
				m.addmatch(als)
				res[v[0].getID()] = [m]
			else:
				continue
		return res		
	except Exception as err:
		print(err)
		sys.exit(0)
# correct support file format
def correctSuppFormat(file):
	handler = open(file,'r')
	result = ''
	for line in handler:
		if len(line.rstrip().split()) == 5 or len(line.rstrip().split()) == 6:
			result += line
			continue
		
		rec = line.rstrip().split('HWI')
		for i in rec:
			if i != '' or i:
				tmp = 'HWI' + i
				result += tmp.rstrip()
				result += '\n'
	handler.close()
	out = open(file,'w')
	out.write(result)
	out.close()	
# only keep those unique alignments with repect to the reads ID
def removeDuplicateSupp(file):
	handler = open(file,'r')
	unique = {}
	for line in handler:
		rec = line.rstrip().split()
		unique[rec[0] + rec[1]] = line
	handler.close()
	handler = open(file,'w')
	for k,v in unique.items():
		handler.write(v)
def parseSuppFileInd(file):
	handler = open(file,'r')
	for line in handler:
		rec = line.rstrip().split()
		al = alignment(rec[0],rec[1],rec[2],rec[3],rec[4],rec[5])
		yield al	
if __name__ == '__main__':
	import glob	
	import os
	dir = '/home/luzhao/THR104_analysis/finalmutation'
	#samples = sampleInfo()
	#sampleName = samples.samples()
	#sampleNames = sampleName['normal'] + sampleName['ET']
	#for e in sampleNames:
	#	if e == 'THR116':
	#file = 'mutation_common_except_normal.supp'
	os.chdir(dir)
	for e in glob.glob('*.supp'):	
		removeDuplicateSupp(os.path.join(dir,e))
			
