import sys
import numpy as np
class seqGenerator(object):
    def __init__(self,pa,pc,pg,pt,pn,prefix,lenOfseq):
        self.__bpgen = bpGenerator(pa,pc,pg,pt,pn)
        self.__filename = prefix + "sequence"
        self.__buf = ""
        self.__lenOfseq = lenOfseq
    def generateStraight(self):
        for i in range(0,self.__lenOfseq):
            self.__buf += self.__bpgen.nextbp()
            f = open(self.__filename,'w')
            f.write(self.__buf)
            f.close()
        return self.__buf
    #need more work
    def generateCircle(self):
        return self.__buf
    #need more work
    def generateBranch(self):
        return self.__buf
class bpGenerator(object):
    def __init__(self,pa,pc,pg,pt,pn):
        self.__prob = [pa,pc,pg,pt,pn]
        if pa < 0 or pc < 0 or pg < 0 or pt < 0 or pn < 0 or sum(self.__prob) != 1:
            print "invalid probability setting"
            sys.exit(0)
        self.__alpha = ['A','C','G','T','N']
    def nextbp(self):
        #print np.random.multinomial(1, self.__prob, size=1).tolist()[0]
        index = np.random.multinomial(1, self.__prob, size=1).tolist()[0].index(1);
        return self.__alpha[index]

