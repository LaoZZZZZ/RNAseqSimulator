#!/usr/bin/python
from parseRefGene import parseRefGene
from parseMutSupp import parseMutSupp
import os
def recalibrate(suppfile,refgene):
	supp = parseMutSupp(suppfile)
	res = []
	for k,v in supp.items():
		for c in v:
			for e in c.getMatch()[0]:
				if not isinstance(e,str):
					d = e.getOrient()
					if refgene.orient() == '-':
						if d == 0:
							e.setOrient(16)
						else:
							e.setOrient(0)
					pos = e.getPos()
					mpos = refgene.getMrnaPos(pos)
					if refgene.orient() == '-':
						e.setPos(mpos-49)
					else:
						e.setPos(mpos)
					res.append(e)
	return res 

if __name__ == '__main__':
	
	inputdir = '/home/luzhao/THR104_analysis/notch2'	
	refgene = parseRefGene('/home/luzhao/THR104_analysis/chromo/refGene.txt')
	mrna = 'NM_001200001'
	refs = refgene[mrna]
	for e in refs:
		if e.rnaname() == mrna+'ex22':
			mrnar = e
	print(mrnar.orient())
	sample = ['A019','A023','N084','THR053','THR039','THR100','THR101','THR102','THR103','THR104','THR116','THR136','THR137']
	for s in sample:
		out = open(os.path.join(inputdir,'notch2_'+s+'_normal_cal.sup'),'w')
		suppfile = os.path.join(inputdir,'notch2_' + s + '_normal.supp')
		result = set(recalibrate(suppfile,mrnar))
		for e in result:
			out.write(e.__str__())
		out.close()					



