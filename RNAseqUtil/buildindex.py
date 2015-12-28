#!/usr/bin/python


def buildindex(inputdir,ref,outdir):
	import os
	f = os.path.join(inputdir,ref)
	out = os.path.join(outdir,ref)
	os.system('bowtie2-build ' + f + '.fa ' + out) 

if __name__ == '__main__':
	import sys
	direct = sys.argv[1]
	ref = sys.argv[2]
	buildindex(direct,ref,direct)	
