#!/usr/bin/python
import os
import sys
from parseMutSupp import parseMutation
from mutation import mutation
from mutation_detect import loadSamByPos
def removeDuplicate(mutfile,suppfile,normfile):
	try:
		mutations = parseMutation(mutfile,suppfile,normfile)
		print(len(mutations))
		uniq = {}
		head = ''
		for mut in mutations:
			key = str(mut.mutationType())+mut.getChromo()+':'+str(mut.getChrMutRange()[0])+'-'+	str(mut.getChrMutRange()[1])
			if key in uniq:
				uniq[key].addSupport(mut.getSupport())
				uniq[key].addNormalSupport(mut.getNormalSupport())
			else:
				uniq[key] = mut	
		basename = os.path.splitext(os.path.basename(mutfile))[0]
		outmutf = open(os.path.join(os.path.split(mutfile)[0],basename+'_unique.csv'),'w')
		outsuppf = open(os.path.join(os.path.split(mutfile)[0],basename+'_unique.supp'),'w')
		outnormsuppf = open(os.path.join(os.path.split(mutfile)[0],basename+'_unique.norm'),'w')
		head = '@mutationType: 1 is insertion 2 is deletiong 3 is snp\n'
		head = head + '@column name\n'
		head = head + '@mutation_type,mrna,mrna_pos,chr,chr_start,chr_end,gene,mutation_length,mutation_seq,uniquesupp,totalsupp,uniquenorm,totalnorm,uniquepercent,totalpercent,supportReadsID,normalReadsID\n'
		outmutf.write(head)
		print(len(uniq))
		for k,e in uniq.items():
			outmutf.write(e.__str__() + '\n')
			for m in e.getSupport():
				outsuppf.write(m.__rec__()+'\n')
			for m in e.getNormalSupport():
				outnormsuppf.write(m.__rec__()+'\n')
		outsuppf.close()
		outmutf.close()
		outnormsuppf.close()
	except Exception as err:
		print(err)
		sys.exit(0)	

if __name__ == '__main__':
	mutfile = '/home/luzhao/THR104_analysis/finalmutation/THR100_small.csv'
	suppfile = '/home/luzhao/THR104_analysis/finalmutation/THR100_small.supp'
	removeDuplicate(mutfile,suppfile)
