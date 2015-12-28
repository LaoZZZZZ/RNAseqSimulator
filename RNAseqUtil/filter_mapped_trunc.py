#! /usr/bin/python

import sys
import os
from util.command.subproc import *
from multiprocessing import Process
def  source_parsing(file):
	dir = ''
	targets=[]
	for line in file:
		line = line.rstrip('\n')
		tuple = line.split(':')
		if tuple[0] =='directory':
			dir=tuple[1]
		elif tuple[0] == 'target':
			targets = targets + tuple[1:]
		else:
			print('undefined name' + tuple[0])
			sys.exit()
	print(dir)
	print(targets)
	return {'dir':dir,'targets':targets}
def run(input,output):
	#command = 'samtools view -F 4 ' +  input +' '  +  ' > ' + output 
	#print(command)	
	#os.popen(command)
	fhandler = open(input,'r')
	outhandler = open(output,'w')
	abnormal = open(output+'_abnormal','w')
	for line in fhandler:
		if line[0] == '@':
			outhandler.write(line)
			continue	
		else:
			record = line.split()
			if record[1] == '4':
				continue
			elif record[1].isdigit() and record[1] != '4':
				outhandler.write(line)
			else:
				abnormal.write(line)
	outhandler.close()
	fhandler.close()
	abnormal.close()	
def filterMapped(inputdir,target,outputdir):
	input = os.path.join(inputdir,target)
	output = os.path.join(outputdir,target+'.sam')
	run(input,output)					
if __name__ == '__main__':
	print('executing')
	task = source_parsing(sys.stdin)
	outputdir = '/home/luzhao/calr_analysis/regine_mapped'
	if not os.path.isdir(outputdir):
		os.system('mkdir '+outputdir)
	dir = task['dir']
	targets=task['targets']
	suffix = ''
	for e in targets:
		tmp = os.path.join(dir,e)
		if os.path.isfile(tmp):
			output = os.path.join(outputdir,e+'.sam')
			Process(target=run,args=(tmp,output)).start()
	print("main process exit")
