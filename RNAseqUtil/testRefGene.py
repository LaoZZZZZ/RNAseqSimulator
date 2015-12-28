#!/usr/bin/python

if __name__ == '__main__':
	import os
	import sys
	from refGeneStruct import refGene
	mrna = refGene('NM_001','chr1',0,30,10,28,'-',3,[3,10,18],[9,15,29],'ABCD')
	seq = 'ATCGGGCAACCCAAAACCCCGGCCGGAAG'
	print(mrna)	
	print(len(seq),mrna.reverseCom(seq))
	print(mrna.getRNASeq(seq),mrna.getCodingRegion(seq))
