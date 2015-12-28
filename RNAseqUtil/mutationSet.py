#!/usr/bin/python
import os
import sys
from parseMutSupp import parseMutation
class mutationSet(object):
	def __init__(self,mutfile):
		if isinstance(mutfile,list):
			muts = parseMutation(mutfile[0],mutfile[1],mutfile[2])
			self.__mutations = {}
			for e in muts:
				key = str(e.mutationType()) + e.getChromo()+':'+str(e.getChrMutRange()[0])+e.mutationSeq()
				if not key in self.__mutations:
					self.__mutations[key] = e
				else:
					self.__mutations[key].addNormalSupport(e.getNormalSupport())
		elif isinstance(mutfile,dict):
			self.__mutations = mutfile
		else:
			raise Exception('invalid variable type %s, only accept filename or dictionary'%(str(type(mutfile))))
	def getMutations(self):
		return self.__mutations	
	def intersect(self,other):
		res = {}
		inter = set(self.__mutations.keys()) & set(other.getMutations().keys())
		for e in inter:
			res[e] = self.__mutations[e]
		return mutationSet(res)	
	def difference(self,other):
		res = {}
		diff = set(self.__mutations.keys()) - set(other.getMutations().keys())
		for e in diff:
			res[e] = self.__mutations[e]
		return mutationSet(res)
	def union(self,other):
		res = {}
		for k,v in self.__mutations.items():
			res[k] = v
		for k,v in other.getMutations().items():
			res[k] = v
		return mutationSet(res)
	def __len__(self):
		return len(self.__mutations)
	def __str__(self):
		head = '@mutationType: 1 is insertion,2 is deletion,3 is snp\n'
		head = head + '@column name\n'
		head = head + '@mutation_type,mrna,mrna_pos,chr,chr_mutation_start,chr_mutation_end,gene,mutation_length,mutation_seq,uniquesupport,totalsupport,uniquenorm,totalnorm,uniquepercent,totalpercent,supportReadsID,normalReadsID\n'
		mut = head
		supp = ''
		norm = ''
		for k,v in self.__mutations.items():
			mut = mut + v.__str__() + '\n'
			for s in v.getSupport():
				supp = supp + s.__rec__() +'\n'
			for s in v.getNormalSupport():
				norm = norm + s.__rec__() + '\n'
		return (mut,supp,norm) 
if __name__ == '__main__':
	dir = '/home/luzhao/THR104_analysis/finalmutation'
	file = 'THR101_unique_annotated.csv'
	mut = mutationSet(os.path.join(dir,file))
	mut2 = mutationSet(os.path.join(dir,'THR102_unique_annotated.csv'))
	print(len(mut.getMutations()),len(mut2.getMutations()))
	print(len(mut.intersect(mut2)),len(mut2.difference(mut)),len(mut.difference(mut2)))
	print(len(mut.union(mut2)))
	muts = [mut,mut2]
	print(len(reduce(mutationSet.difference,muts)))
