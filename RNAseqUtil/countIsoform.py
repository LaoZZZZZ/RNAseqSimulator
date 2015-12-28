#!/usr/bin/python
from parseRefGene import parseRefGeneByGene
def countIsoform(refgenefile):
	
	refgene = parseRefGeneByGene(refgenefile)
	result = {}
	for k,v in refgene.items():
		if len(v) in result:
			result[len(v)] += 1
		else:
			result[len(v)] = 1
	return result

if __name__ == '__main__':
	file = '/home/luzhao/THR104_analysis/chromo/refGene_new.txt'
	res = countIsoform(file)
	outfile = open('/home/luzhao/isoform_stat.csv','w')
	for k,v in res.items():
		outfile.write('%s,%s\n'%(k,v))
	outfile.close()
