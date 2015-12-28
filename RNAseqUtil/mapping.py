#!/usr/bin/python

import os
import sys
import subprocess 

from extract_mapped import source_parsing

def pairmapping(inputdir,fq1,fq2,genome,outputdir,prefix,k =20,gapped = 20, mis = 1):
	try:
		program = 'bowtie2'
		options = '-N ' + str(mis) + '  --no-hd  -p 8 --gbar ' + str(gapped) + ' -k ' + str(k)  
		inputfq1 = os.path.join(inputdir,fq1)
		inputfq2 = os.path.join(inputdir,fq2)
		if not os.path.isfile(inputfq1) or not os.path.isfile(inputfq2):
			raise Exception('%s and %s is not reads file'%(inputfq1,inputfq2))
		#if not unmappout:
		#	unmappout = outputdir
		outputfile = os.path.join(outputdir,prefix+'.sam')
		#if repUn:
		#	unmapped = os.path.join(unmappout,prefix+'_unmapped.fq')
	      	#	command = program + ' ' + options + ' ' + genome + ' -q  ' + ' -1 ' +  inputfq1 + ' -2 ' + inputfq2 + ' -S ' +  outputfile  + ' --un ' + unmapped
		#else:
	      	command = program + ' ' + options + ' ' + genome + ' -q  ' +  ' -1 ' +  inputfq1 + ' -2 ' + inputfq2 + ' -S ' +  outputfile 
		print(command) 	
		pipe = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
		pipe.communicate()
		if pipe.returncode:
			raise Exception('Something wrong when doing the alighment using' + program )	 
	except Exception as error:
		print(error)
		sys.exit(0)
def mapping(inputdir,target,genome,outputdir,outfile = '',unmappout = '',repUn=True,k=10):
	try:
		program = 'bowtie2'
		options = '-N 1 --no-unal --no-hd  -p 2 ' + ' -k ' + str(k)
		inputfq = os.path.join(inputdir,target+'.fq')
		if not os.path.isfile(inputfq):
			raise Exception('%s is not reads file'%(inputfq))
		if not unmappout:
			unmappout = outputdir
		if not outfile:
			outfile = target
		outputfile = os.path.join(outputdir,outfile+'.sam')
		if repUn:
			unmapped = os.path.join(unmappout,target+'_unmapped.fq')
	      		command = program + ' ' + options + ' ' + genome + ' -U  ' + inputfq + ' -S ' +  outputfile  + ' --un ' + unmapped
		else:
	      		command = program + ' ' + options + ' ' + genome + ' -U  ' + inputfq + ' -S ' +  outputfile 
		print(command) 	
		pipe = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
		pipe.communicate()
		if pipe.returncode:
			raise Exception('Something wrong when doing the alighment using' + program )	 
	except Exception as error:
		print(error)
		sys.exit(0)
def gzmapping(inputdir,target,genome,outputdir,p = 8,k = 20,outfile = '',unmappout = '',repUn=True):
	try:
		program = 'bowtie2'
		print(genome,k)
		options = '-N 1 --no-unal --no-hd  -p ' + str(p)   + '  -k ' + str(k)
		inputgz = ''		
		for f in os.listdir(inputdir):
			if f.split('.')[-1] == 'gz' and  'sequence' in f:
				inputgz = inputgz +  os.path.join(inputdir,f) + ','
		#inputfq = os.path.join(inputdir,target+'.fq')
		if not unmappout:
			unmappout = outputdir
		if not outfile:
			outfile = target
		outputfile = os.path.join(outputdir,outfile+'.sam')
		if repUn:
			unmapped = os.path.join(unmappout,target+'_unmapped.fq')
	      		command = program + ' ' + options + ' ' + genome + ' -U  ' + inputgz + ' -S ' +  outputfile  + ' --un ' + unmapped
		else:
	      		command = program + ' ' + options + ' ' + genome + ' -U  ' + inputgz + ' -S ' +  outputfile 
		print(command) 	
		pipe = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
		pipe.communicate()
		if pipe.returncode:
			raise Exception('Something wrong when doing the alighment using' + program )	 
	except Exception as error:
		print(error)
		sys.exit(0)
if __name__ == '__main__':
	print('mapping')
	import os
	import sys
	data = source_parsing(sys.stdin)
	direct = data['dir']
	genome = data['genome']
	targets = data['targets']
	print(direct,genome,targets)
	output = '/home/luzhao/multipleHits/'
	from multiprocessing import Process
	procs = []
	import time
	#time.sleep(3600)
	for target in targets:
		p = Process(target = gzmapping,args = (os.path.join(direct,target),target,genome,output))
		p.start()
		procs.append(p)
	for p in procs:
		p.join()	
