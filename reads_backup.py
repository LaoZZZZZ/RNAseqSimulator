import math
import sys
from sets import Set
import numpy as np
import pdb
# class holds single reads    
class read(object):
	def __init__(self,id,seq,qua):
		self.__id = id
		self.__seq = seq
		self.__qua = qua
		if len(self.__seq) != len(self.__qua):
			print "the length of quality is not equal to the length of read sequence!"
			sys.exit(0)
	def __len__(self):
		return len(self.__seq)
	def seq(self):
		return self.__seq
	def name(self):
		return self.__id
	def qual(self):
		return self.__qua
	def empty(self):
		return not self.__seq
	def __repr__(self):
		return '%s\n%s\n+\n%s'%(self.__id,self.__seq,self.__qua)
	def __str__(self):
		return self.__repr__()
 		
    # find the reverse complementary given a sub seq
    def revComplementary(self,seq):
        rc = ""
        for c in reversed(seq):
            if c in ('A','a'):
                rc += 'T'
            elif c in ('C','c'):
                rc += 'G'
            elif c in ('G','g'):
                rc += 'C'
            elif c in ('T','t'):
                rc += 'A'
            elif c in ('N','n'):
                rc += 'N'
            else:
                print "invalid sequence!"
	return rc

class rdGenerator(object):
    def __init__(self,seq,coverage,lenOfread,lenOfkmer,prefix,form = "fq"):
        self.__seq = seq
        self.__cov = coverage
        self.__lenOfread = lenOfread
        self.__prefix = prefix
        #self.__reads = []
        self.__format = form
        self.__total = int(math.ceil(coverage * len(self.__seq)/lenOfread))
        self.__kmer = {}
        self.__lenOfk = lenOfkmer
        self.__generatedReads = Set()
        self.__uniqKmer = []  #sorted unique pattern exist in the whole sequence
        self.__kdict=[] #  # keep track of each kmer at the ith position of original sequence
                           # the index of this array is corresponding to the position of original sequence
                           # the value at each index corresponding to the index of pattern in the self.__uniqKmer
        self.__krevdict=[] # keep track reverse complementary of each kmer at the ith position of original sequence
                           # the index of this array is corresponding to the position of original sequence
                           # the value at each index corresponding to the index of pattern in the self.__uniqKmer
        self.__kstat = {}  # keep track of the frequenc of each kmer that shows in the reads
    # interface
    def generate(self):
        self.__findUniqueKmer()
        self.__output()

    # private functions
    def __findUniqueKmer(self):
        tmp = Set()
        for i in range(len(self.__seq) - self.__lenOfk + 1):
            if self.__seq[i: (i + self.__lenOfk)] not in tmp:
                tmp.add(self.__seq[i: (i + self.__lenOfk)])
                tmp.add(self.__revComplementary(self.__seq[i: (i + self.__lenOfk)]))
        self.__uniqKmer = list(tmp)
        self.__uniqKmer.sort()
        for i in range(len(self.__uniqKmer)):
            self.__kstat[i] = 0
        #print max(tmp),len(tmp),min(tmp)
        tmp = Set()
        #self.__uniqKmer = Set([])
        for i in range(len(self.__seq) - self.__lenOfk + 1):
            self.__kdict.append(self.__uniqKmer.index(self.__seq[i: (i + self.__lenOfk)]))
            self.__krevdict.append(self.__uniqKmer.index(self.__revComplementary(self.__seq[i: (i + self.__lenOfk)])))
        
    def __generateRead(self):
        index = np.random.uniform(0,len(self.__seq) -1,1)
        index = int(index)
        attempt = 0
        index = self.__calIndex(index)
        # guarantee that there is no repeated reads
        while index in self.__generatedReads:
            index = np.random.uniform(0,len(self.__seq) -1,1)
            index = int(index)
            attempt += 1
            if attempt > 50:
                return read(0,"","")
            index = self.__calIndex(index)
       # pdb.set_trace()
        self.__generatedReads.add(index)
        return read(index,self.__seq[index:index+self.__lenOfread ],self.__qual())
        
        
            
    def __writeReads(self,r,fr):
        if self.__format == "fq":
            fr.write('@{0}\n{1}\n+\n{2}\n'.format(r.name(),r.seq(),r.qual()))
        elif self.__format == "fa":
            fr.write('>{0}\n{1}\n'.format(r.name(),r.seq()))
        else:
            print "only support fq and fa reads file"
            sys.exit(0)
        #print "writing reads"
    def __updateKmers(self,r):
        for i in range(self.__lenOfread - self.__lenOfk + 1):
            if not self.__kstat.has_key(self.__kdict[i + r.name()]):
                print "something wrong when generating the kmer"
                sys.exit(0)
            #pdb.set_trace()
            self.__kstat[self.__kdict[i + r.name()]] += 1
            self.__kstat[self.__krevdict[i + r.name()]] += 1
            
    def __output(self):
        fr = open(self.__prefix + "reads." + self.__format,"w")
        self.__basicCover(fr)
        for i in range(1,self.__total):
            r = self.__generateRead()
            if r.empty():
                break
            self.__writeReads(r,fr)
            self.__updateKmers(r)
        fr.close()

        if len(self.__kstat) != len(self.__uniqKmer):
            print "something wrong when generating the kmers\n"
            sys.exit(0)
        self.__writeKmers()
    def __basicCover(self,fr):
        for i in range(len(self.__seq) - self.__lenOfread + 1):
            r = read(i,self.__seq[i:i+self.__lenOfread],self.__qual())
            self.__writeReads(r,fr)
            self.__updateKmers(r)
            self.__generatedReads.add(i)
            if min(i + self.__lenOfk - self.__lenOfk,len(self.__seq) - self.__lenOfread) < i:
                break                
            i = min(i + self.__lenOfk - self.__lenOfk,len(self.__seq) - self.__lenOfread) 
    # write kmer to the file
    def __writeKmers(self):
       fk = open(self.__prefix + "kmer","w")
       for k, v in self.__kstat.iteritems():
           fk.write('{0}\t\t{1}\n'.format(self.__uniqKmer[k],v))
       fk.close()
    def __qual(self):
        s = ""
        for i in range(0,self.__lenOfread):
            s += chr(int(np.random.uniform(33,126,1)))
        return s
    # decide the start point given the random index
    def __calIndex(self,index):
        if index < self.__lenOfread - 1 and len(self.__seq) - index < self.__lenOfread:
            if self.__lenOfread > 2 * index:
                index = 0
            elif self.__lenOfread > 2*(len(self.__seq) - index):
                index = len(self.__seq) - self.__lenOfread
            else:
                index = index - int(self.__lenOfread/2)
        elif index >= self.__lenOfread - 1:
            index = index - self.__lenOfread + 1
        elif len(self.__seq) - index + 1>= self.__lenOfread:
            index = index
        else:
            print "invalid sequence length({0}) vs length of read({1})".format(len(self.__seq),self.__lenOfread)
            sys.exit(0)
        return index
    # find the reverse complementary given a sub seq
    def __revComplementary(self,seq):
        rc = ""
        for c in reversed(seq):
            if c in ('A','a'):
                rc += 'T'
            elif c in ('C','c'):
                rc += 'G'
            elif c in ('G','g'):
                rc += 'C'
            elif c in ('T','t'):
                rc += 'A'
            elif c in ('N','n'):
                rc += 'N'
            else:
                print "invalid sequence!"
                sys.exit(0)
        Return rc
