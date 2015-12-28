#!/usr/bin/python
def parse_input(file):
	res = {}
	for line in file:
		item = line.rstrip('\n').split(':')
		if item[0] == 'genome':
			res['genome']=item[1]
		elif item[0] == 'directory':
			res['directory'] = item[1]
		elif item[0] == 'target':
			res['target'] = item[1].split(':')
		elif item[0] == 'options':
			res['options'] = item[1]
		elif item[0] == 'output':
			res['output'] = item[1]
		else:
			print('invalid items \n')
			sys.exit(0)
	return res

def tophat(dir,genome,option,target,output):
	workdir = os.path.join(dir,target)
	files = [] 
	for f in os.listdir(workdir):
		absfile = os.path.join(workdir,f)
		if os.path.isfile(absfile) and f.endswith('.gz') and not 'sequence' in f:
			files.append(absfile)
	inputs = ','.join(files)  
	command = 'tophat ' + option + ' ' + genome + ' ' + inputs 
	print(command)
	pipe = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
	print('start running')
	pipe.communicate()
	attachment = ['tophat_log']
	subject = 'Tophat finished'
	body = 'Tophat has finished the alignment'
	if pipe.returncode:
		subject = 'Error of tophat'
		body = 'Something wrong with the alignment using tophat!\n'
		raise Exception('Something wroing when doing the alignment using tophat \n' + command + '\n')
	email2me(subject,body,attachment)
				
if __name__== '__main__':
	import os
	import sys
	import subprocess
	from multiprocessing import Process
	from myemail import email2me 
	task = parse_input(sys.stdin)
	print(task)
	inputdir = task['directory']
	targets = task['target']
	output = task['output']
	options = task['options']
	genome = task['genome']
	for target in targets:
		Process(target=tophat,args=(inputdir,genome,options,target,output)).start()
	print('Main process exit \n')	
