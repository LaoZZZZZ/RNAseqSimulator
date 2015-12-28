#!/usr/bin/python

from fragmentLength import fragmentLength as fragConfig
from mrnaDetermineGenerator import mrnaDetermineGenerator
from possionProcess import PossionProcess as pp
from fragment import Fragment
from fragmentFilter import filterFragment as filterEngine
from amplification import amplifier


# given mrnalibrary and possion process coming rate,the valid length interval for fragment
# generate all valid fragment 
class fragmentGenerator(object):
	def __init__(self,mrnaGen,lambd,flower,fupper):
		self.mrnaGen = mrnaGen
		self.lambd = lambd
		self.filterEng = filterEngine(flower,fupper)
		self.generateFragment()
	# fragment all full length transcripts by possion process 
	def generateFragment(self):
		print('In the process of fragmentation')
		self.fragments = []
		for m in self.mrnaGen:
			tmp = []
			cutPoints = pp.cutPositions(self.lambd,len(m[1]))
			start = 0
			for c in cutPoints:
				tmp.append(Fragment(m[0],start,c-start))
				start = c
			tmp = self.filterEng.filterFragments(tmp)
			self.fragments += tmp
		# once transcripts fragmented, all fragment will be amplified by an amplifier
		#self.amplifier.amplify(self.fragments)
		self.cur = 0
		print('totally have %s fragments'%(len(self.fragments)))
		print('Fragmentation finished!')				
	def __iter__(self):
		return self
	def next(self):
		if self.cur < len(self.fragments):
			tmp = self.cur
			self.cur+=1
			return self.fragments[tmp]
		else:
			raise StopIteration()
	def __len__(self):
		return len(self.fragments)
if __name__ == '__main__':
	mrnafile = 'mrna.txt'
	lambd = 10 
	total = 1
	generator = fragmentGenerator(mrnaDetermineGenerator(mrnafile,10),lambd,1,100)
	for f in generator:
		print(f)
        
