#!/usr/bin/python
import os
from mutationSet import mutationSet
from parseMutSupp import correctSuppFormat,removeDuplicateSupp
from sampleInfo import sampleInfo
def find_common(prefix,common_set,outprefix):
	mut = prefix + '.csv'
	supp = prefix + '.supp'
	norm = prefix + '.norm'
	mutations = mutationSet([mut,supp,norm]).intersect(common_set)

	mut,supp,norm = mutations.__str__()	
	outmut = open(outprefix+'.csv','w')
	outmut.write(mut)
	outmut.close()
	outsupp = open(outprefix+'.supp','w')
	outsupp.write(supp)
	outsupp.close()
	correctSuppFormat(outprefix+'.supp')
	removeDuplicateSupp(outprefix+'.supp')
	outnorm = open(outprefix+'.norm','w')
	outnorm.write(norm)
	outnorm.close()

if __name__ == '__main__':
	samples = sampleInfo()
	ets = samples.norm()
	directory = '/home/luzhao/THR104_analysis/finalmutation'
	outdir = os.path.join(directory,'possible_mutations')
	commfiles = [os.path.join(directory,'mutation_common_except_normal.csv'),os.path.join(directory,'mutation_common_except_normal.supp'),os.path.join(directory,'mutation_common_except_normal.norm')]
	comm = mutationSet(commfiles)
	for e in ets:
		prefix = directory + '/' + e + '_small_unique'
		outprefix = outdir + '/' + e + '_small'	
		find_common(prefix,comm,outprefix)	
