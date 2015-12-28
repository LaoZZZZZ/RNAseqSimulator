#!/usr/bin/python
# split a read into two read at the position pos , pos start at 0
# split result  [0,pos],[pos+1,end]
import os
import sys
def split(pos,rd):
	from reads import read
	if pos < 0 :
		print('the position ({0}) is negative'.format(pos))
	elif len(rd) < pos - 1:
		print('the length({0}) of read is smaller than the split position({1}'.format(len(rd),pos))
		return rd
	else:
		rd1 = read(rd.getid(),rd.getseq()[0:pos+1],rd.getqual()[0:pos+1])
		rd2 = read(rd.getid(),rd.getseq()[pos+1:],rd.getqual()[pos+1:])
		return [rd1,rd2]
def getFrontEnd(sublen,rd):
	from reads import read
	if sublen < 0 :
		print('the length({0}) is negative'.format(sublen))
	elif len(rd) < sublen:
		print('the length({0}) of read is smaller than the trancated length({1}'.format(len(rd),sublen))
		return rd
	else:
		length = len(rd)
		rd1 = read(rd.getid(),rd.getseq()[0:sublen],rd.getqual()[0:sublen])
		rd2 = read(rd.getid(),rd.getseq()[length-sublen:],rd.getqual()[length - sublen:])
		return [rd1,rd2]
def loadOneRead(fhandler):
	try:
		from reads import read
		rd = []
		total = 0
		for line in fhandler:
			total = total + 1
			rd.append(line.rstrip('\n'))
			if total == 4:
				if rd[2] != '+':
					raise Exception('invalid fq file format!\n')
				else:
					break
		if rd:
			return read(rd[0],rd[1],rd[3])
		else:
			return None
	except Exception as err:
		print(err)
		import sys
		sys.exit(0)
def writeOneRead(fhandler,rd):
	try:
		fhandler.write(rd.__str__())
		fhandler.write('\n')
	except Exception as err:
		print(err)
		import sys
		sys.exit(0)
def writeReads(out,rds):
	try:
		from reads import read
		for rd in rds:
			writeOneRead(out,rd)
	except Exception as err:
		print(err)
		import sys
		sys.exit(0)
def readIterate(file):
	fhandler = open(file,'r')
	rd = loadOneRead(fhandler)
	while rd:
		yield rd
		rd = loadOneRead(fhandler)	

def splitreads(file,output,sublen):
        from reads import read
        out = open(output,'w')
        for rd in readIterate(file):
                rds = getFrontEnd(sublen,rd)
                writeReads(out,rds)
        out.close()
        fhandler.close()
def splitReadsPair(file,outdir,prefix,sublen):
        from reads import read
        fhandler = open(file,'r')
        head = open(os.path.join(outdir,prefix+'_1.fq'),'w')                              
        tail = open(os.path.join(outdir,prefix+'_2.fq'),'w')
        for rd in readIterate(file):
                rds = getFrontEnd(sublen,rd)
                rd1 = read(rds[0].getid().split('/')[0] + '/1',rds[0].getseq(),rds[0].getqual())
                rd2 = read(rds[1].getid().split('/')[0] + '/2',rds[1].getseq(),rds[1].getqual())
                writeOneRead(head,rd1)
                writeOneRead(tail,rd2)
        head.close()
        tail.close()

if __name__ == '__main__':
	
	from reads import read
	print('this is reads functions package')
	dirc = '/home/luzhao/THR104_analysis/unmapped'
	ofq = os.path.join(dirc,'THR104.fq')
	splitReadsPair(ofq,dirc,'THR104',20)
