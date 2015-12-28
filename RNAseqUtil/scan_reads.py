#! /usr/bin/python
import gzip
import os
import sys
from multiprocessing import Process,Lock

from extract_mapped import source_parsing
from reads import read
from read_func import split, getFrontEnd


def trim(id):
	tmp = id.split('/')[0]
	if tmp[0] == '@':
		tmp = tmp[1:]
	return tmp
def loadSam(samdir,target,start = 0,end = float('Inf')):
	from sets import Set
	samfile = os.path.join(samdir,target+'.sam')
	print(samfile)
	dict = Set()
	if not os.path.isfile(samfile):
		print('{0} is not a sam file'.format(samfile))
		sys.exit()
	fhandle = open(samfile,'r')
	for line in fhandle:
		if line[0] != '@':
			rec = line.rstrip('\n').split()
			startpos = int(rec[3])
			if startpos >= start and startpos <= end:		
				dict.add(trim(rec[0]))
	return dict
def singleSample(gzdir,samdir,target,output,sublen,dataset=0):
	workdir = os.path.join(gzdir,target)
	dict = loadSam(samdir,target)		
	outputfile = os.path.join(output,target+'.unmapped')
	outhandler = open(outputfile,'w')
	totalreads = 0.0 
	totalunmapped = 0.0
	totalmapped = 0
	originalunmapped = open(os.path.join(output,target+'.fq'),'w')
	counted_gz = []
	sequence_gz = []
	for f in os.listdir(workdir):
		absfile = os.path.join(workdir,f)
		if os.path.isfile(absfile) and absfile.endswith('.gz') and not 'sequence' in f:
			counted_gz.append(absfile)
		elif os.path.isfile(absfile) and absfile.endswith('.gz') and 'sequence' in f:
			sequence_gz.append(absfile)
	candidate = []	
#   choose the counted dataset, 0 means prefer counted reads files, 1 means prefer 'sequence' reads file	
	if dataset == 1:
		candidate = counted_gz 
	else:
		candidate = sequence_gz 
	for e in candidate:
		print(e)
		file = gzip.open(e,'r')
		total = 0
		rd = []	
		for line in file:
			total = total + 1
			rd.append(line.rstrip('\n'))
			if total %  4 == 0:
				rd[0] = rd[0].split()[0]
				if not '+' in rd[2]:
					print('incorrect fq gz format')
					sys.exit(0)
				rdbody = read(rd[0],rd[1],rd[3])
				totalreads = totalreads + 1
				if not trim(rd[0]) in dict: 
					totalunmapped = totalunmapped + 1
					originalunmapped.write(rdbody.__str__())
					originalunmapped.write('\n')
					#print(rdbody)
					splitrd = getFrontEnd(sublen,rdbody)
					outhandler.write(splitrd[0].__str__())
					outhandler.write('\n')
					outhandler.write(splitrd[1].__str__())
					outhandler.write('\n')
				else:
					totalmapped = totalmapped + 1	
				rd = []
	if totalreads > 0:
		stat = open(os.path.join(output,target+'.stat'),'w')
		stat.write('Total : %s \n Mapped reads: %s \n Unmapped reads: %s  \n Fraction of mapped reads: %s \n Fraction of  unmapped: %s' % (totalreads,totalmapped,totalunmapped,totalmapped/totalreads,totalunmapped/totalreads))
		stat.close()	
	originalunmapped.close()
	outhandler.close()
def hello():
	print(os.getpid())
if __name__ == '__main__':
	task = source_parsing(sys.stdin)
	print(task)
	samdir = '/home/luzhao/calr_analysis/samfile'
	output = '/home/luzhao/calr_analysis/unmapped'
	if not os.path.isdir(output):
		os.system('mkdir ' + output)		
	for target in task['targets']:
		print(target)
		#Process(target=hello,args=()).start()
		#Process(target=singleSample,args=(task['dir'],samdir,target,output,20)).start()
		singleSample(task['dir'],samdir,target,output,20,0)
		#dict = loadSam(samdir,target)
		#print(len(dict))
	print('finished filtering')
