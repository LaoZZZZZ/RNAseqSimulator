import hashlib
import binascii
class read:
        hashlen = 30
        Nstrat = 'A'
        dictmap = {'A':0,'a':0,'C':3,'c':3,'G':2,'g':2,'T':1,'t':1}
        def __init__(self,ID,seq,quality):
                try:
                        self.id = ID
                        self.seq=seq
                        self.qual=quality
                        if quality and len(seq)!=len(quality):
                                raise Exception('unequal length of seqence and quality')
                        self.len = len(seq)
                except Exception as err:
                        print(err)
        def getid(self):
        	return self.id
        def getseq(self):
        	return self.seq
        def getqual(self):
        	return self.qual
        def bitwise(rd):
        	fw = 0
        	itr = 0
        	seq = rd.getseq()
        	while itr < read.hashlen:
        		c = seq[itr]
        		fw = fw << 2
        		if c in read.dictmap:
        			fw |= read.dictmap[c]
        		elif c in ['n','N']:
        			fw |= read.dictmap[read.Nstrat]
        		else:
        			raise Exception('There are some other character other than A,C,G,T,N which is ' +c)
        		itr = itr + 1
        	return fw
        def __repr__(self):
        	return '{0}\n{1}\n+\n{2}'.format(self.id, self.seq ,self.qual)	
        def __str__(self):
        	return self.__repr__()
        def __len__(self):
        	return self.len
        def __hash__(self):
                return int(hashlib.md5(binascii.a2b_qp(self.id)).hexdigest(),16)
        def __eq__(self,other):
                return self.seq == other.getseq()
