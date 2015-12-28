#!/usr/bin/python
import os
import sys
from mutationSet import mutationSet
from annotation import countMapped
from mutation_detect import trim,loadSamImp
from extract_mapped import extractRegionMapping,extractSamRegion
def fqfileCount(file):
	h = open(file,'r')
	total = 0
	for line in h:
		total = total + 1
	return total / 4
def  smallSamplesCompile(inputdir,inputfqdir,samples,prefix):
	smalls = []
	out = open(os.path.join(inputdir,prefix + '.csv'),'w')
	for  e in samples:
		small = mutationSet([os.path.join(inputdir,e + '_small_unique.csv'),os.path.join(inputdir,e + '_small_unique.supp'),os.path.join(inputdir,e + '_small_unique.norm')])
		smalls.append(small)
	head = '@' +'sample' + ',' +  'small' +  ',' + 'small_unique' + ','  + 'small_common'   '\n'
	out.write(head)
	small_common = reduce(mutationSet.intersect,smalls)
	for i in range(len(smalls)):
		small_unique = smalls[i].difference(reduce(mutationSet.union,[e for e in smalls if e != smalls[i]]))
		out.write('%s,%s,%s,%s\n'%(samples[i],len(smalls[i]),len(small_unique),len(small_common)))
		out.close()
def  samplesCompile(inputdir,inputfqdir,samples,prefix):
	smalls = []
	larges = []
	overalls = []
	#out = open(os.path.join(inputdir,prefix + '.csv'),'w')
	for  e in samples:
		small = mutationSet([os.path.join(inputdir,e + '_small_unique.csv'),os.path.join(inputdir,e + '_small_unique.supp'),os.path.join(inputdir,e + '_small_unique.norm')])
		#large =	mutationSet(os.path.join(inputdir,e + '_large_unique_annotated.csv'))
		#overall = mutationSet(os.path.join(inputdir,e + '_unique_annotated.csv'))
		smalls.append(small)
		#larges.append(large)
		#overalls.append(overall)
	head = '@' +'sample' + ',' +  'small' + ',' + 'medium/large' + ',' + 'overall' + ',' + 'small_unique' + ',' + 'large_unique' +',' +  'overall_unique' + ',' + 'small_common' + ',' + 'large_common' + ',' + 'overall_common' + ',' + 'unmapped_reads' +  '\n'
	#out.write(head)
	#overall_common = reduce(mutationSet.intersect,overalls)
	#large_common = reduce(mutationSet.intersect,larges)
	small_common = reduce(mutationSet.intersect,smalls)
	#for i in range(len(smalls)):
	#	unique = overalls[i].difference(reduce(mutationSet.union,[e for e in overalls if e != overalls[i]]))
	#	unique_h = open(os.path.join(inputdir,'unique_' + samples[i]+'.csv'),'w')
		#print(unique.__str__())
	#	unique_h.write(unique.__str__())
	#	small_unique = smalls[i].difference(reduce(mutationSet.union,[e for e in smalls if e != smalls[i]]))
	#	large_unique = larges[i].difference(reduce(mutationSet.union,[ e for e in larges if e !=larges[i]]))
	#	out.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(samples[i],len(smalls[i]),len(larges[i]),len(overalls[i]),len(small_unique),len(large_unique),len(unique),len(small_common),len(large_common),len(overall_common),fqfileCount(os.path.join(inputfqdir,samples[i] + '.fq'))))
	#	unique_h.close()
	#print(len(overall_common))
	#commonh = open(os.path.join(inputdir,prefix + '_common.csv'),'w')
	#commonh.write(overall_common.__str__())
	#commonh.close()
	#large_commonh = open(os.path.join(inputdir,prefix + '_large_common.csv'),'w')
	small_commonh = open(os.path.join(inputdir,prefix + '_small_common.csv'),'w')
	small_common_supp = open(os.path.join(inputdir,prefix + '_small_common.supp'),'w')
	small_common_norm = open(os.path.join(inputdir,prefix + '_small_common.norm'),'w')
	#large_commonh.write(large_common.__str__())
	mut,supp,norm = small_common.__str__()
	small_commonh.write(mut)
	#large_commonh.close()
	small_common_supp.write(supp)
	small_common_norm.write(norm)
	small_commonh.close()
	small_common_supp.close()
	small_common_norm.close()

def exceptNorm(inputdir,normals,ets):
	et_mut = mutationSet([os.path.join(inputdir,ets[0]),os.path.join(inputdir,ets[1]),os.path.join(inputdir,ets[2])])
	normal_muts = []
	for e in normals:
		mutfile = os.path.join(inputdir,e+'_small_unique.csv')
		suppfile = os.path.join(inputdir,e+'_small_unique.supp')
		normfile = os.path.join(inputdir,e+'_small_unique.norm')
		normal_muts.append(mutationSet([mutfile,suppfile,normfile]))
	union_normal = reduce(mutationSet.union,normal_muts)
	diff = et_mut.difference(union_normal)
	return diff	
def combineMutations(inputdir,common,samples,outdir,prefix):
	from parseMutSupp import parseMutSupp
	comm = mutationSet(os.path.join(inputdir,common))
	out = open(os.path.join(outdir,prefix+'.csv'),'w')
	head = '@mutationType: 1 is insertion, 2 is deletion  3 is snp\n'
	head = head + '@column name\n'
	head = head + '@sample,mutation_type,mrna,mrna_pos,chr,chr_start,chr_end,gene,mutation_length,mutation_seq,uniquesupport,totalsupport,uniquenorm,totalnorm,uniquepercent,totalpercent,readsID\n'
	out.write(head)
	samplemutation = {}
	for e in samples:
		tmp = mutationSet(os.path.join(inputdir,e+'_unique_annotated.csv'))
		samplemutation[e] = tmp.intersect(comm).getMutations()	
		assert(len(samplemutation[e]) == len(comm))
	for k,v in comm.getMutations().items():
		for s in samples:
			if 'NOTCH2' in  samplemutation[s][k] and '120612005' in samplemutation[s][k]:
				supp = loadSamImp(os.path.join('/home/luzhao/THR104_analysis/mapped',s+'.sam'))
				supph = open(os.path.join(outdir,prefix+'_'+s+'.supp'),'w')
				out.write(s + ',' + samplemutation[s][k])
				ids = samplemutation[s][k].split(',')[15].split()
				for i in ids:
					als = set(supp[trim(i)])
					for a in als:
						supph.write(a.__str__())
				supph.close()
						
					
#		out.write('\n'*4)
	out.close()
def countNormal(inputdir,normals,regions):
	outfile = '/home/luzhao/THR104_analysis/tmp.sam'
	res = {}
	for r in regions:
		for e in normals:
			if r[0] in res:
				res[r[0]].append(e + ',' + r[0] + ',' +  '0'+',' +'0' + ',' + ','.join(map(str,countMapped(os.path.join(inputdir,e+'/'+e+'_accepted_hits.bam'),r[1],outfile))) + '\n')

			else:
				res[r[0]] = [e + ',' + r[0] +','+ '0'+ ',' + '0' + ',' + ','.join(map(str,countMapped(os.path.join(inputdir,e+'/'+e+'_accepted_hits.bam'),r[1],outfile))) + '\n']
	return res
def extractNormal(inputdir,normals,regions,outdir,prefix):
	for r in regions:
		for e in normals:
			outh = open(os.path.join(outdir,prefix+'_' + e + '_normal.supp'),'w')
			extractRegionMapping(os.path.join(inputdir,e+'/'+e+'_accepted_hits.bam'),r[1],'tmp.sam')
			supp = loadSamImp('tmp.sam')
			for k,v in supp.items():
				for  a in v:
					outh.write(a.__str__())
if __name__ == '__main__':
	dir = '/home/luzhao/THR104_analysis/finalmutation'
	outdir = '/home/luzhao/THR104_analysis/notch2'
	inputdir = '/home/stc/Documents/GenomicData/RNAseq/WBahou'
	fqdir = '/home/luzhao/THR104_analysis/unmapped'
	samples = ['THR100','THR101','THR102','THR104','THR116','THR136','THR137']
	#smallSamplesCompile(inputdir,inputfqdir,normals,prefix)
	normals = ['A019','A023','N084','THR039','THR053']
	#ets = 'mutation_common_except_normal.csv'
	#muts = mutationSet(os.path.join(dir,ets))
	#region = []
	#for k,v in muts.getMutations().items():
	#	if '120612005' in v:
	#		rec = v.split(',')
	#		region.append((','.join(rec[0:9]),rec[3]+':'+rec[4]+'-'+rec[5]))
	#extractNormal(inputdir,samples,region,outdir,'notch2')
	#out = open(os.path.join(dir,'normal.csv'),'w')
		
	#for k,v in res.items():
		#for e in v:
		#	out.write(e)
		#out.write('\n'*3)
	samplesCompile(dir,fqdir,normals,'normal_mutation_compile')
	#exceptNorm(inputdir,normals,ets)
	mut,supp,norm = (exceptNorm(dir,samples,['normal_mutation_compile_small_common.csv','normal_mutation_compile_small_common.supp','normal_mutation_compile_small_common.norm'])).__str__()
	outh = open(os.path.join(dir,'normal_mutation_common_except_et.csv'),'w')
	outh.write(mut)
	outh = open(os.path.join(dir,'normal_mutation_common_except_et.supp'),'w')
	outh.write(supp)
	outh = open(os.path.join(dir,'normal_mutation_common_except_et.norm'),'w')
	outh.write(norm)
	outh.close()
	#combineMutations(dir,'mutation_common_except_normal.csv',samples,outdir,'notch2')
			
