import os
import sys
import getopt
import pdb
from simulator import run 
print("Current path is ")
p = os.getcwd()
print p
sys.path.append(p)
def usage():
    print "\tthere are TEN options:"
    print " \t once you specify the probability of any one of A,G,C T and N, you have \
to specify other three, or you can leave it as default(0.25 for each, 0 for N)"
    print "\t A/a(default 0.25,optional): the probability of occurrence of A"
    print "\t C/c(default 0.25,optional): the probability of occurrence of C"
    print "\t G/g(default 0.25,optional): the probability of occurrence of G"
    print "\t T/t(default 0.25,optional): the probability of occurrence of T"
    print "\t N(default 0.0,optional):   the probability of occurrence of N"

    print "\t k(default 22,optional): The length of kmer"
    print "\t l(default 50,optional): The length of reads"
    print "\t v(default 5,optional): coverage of each base"
    print "\t s(required, should be larger or equal len of read): the total length of whole sequence"
    print "\t o(default ./output_,optional): the prefix of output files"
def driver():
    pa = 0.25
    pc = 0.25
    pg = 0.25
    pt = 0.25
    pn = 0
    coverage = 5; # number of reads
    lenOfreads = 50; # length of reads
    lenOfkmer  = 22; # length of kmer
    prefix = os.getcwd() + "/output_"
    lenOfseq = 0 # length of whole sequence
    #pdb.set_trace()
    try:
        opt,args = getopt.getopt(sys.argv[1:],"N:T:t:G:g:C:c:A:a:k:l:v:s:o:h",["help"])
    except getopt.error,msg:
        print msg
        usage()
        print "for help use --help"
        sys.exit(2)
    for o,a in opt:
        if o in ("-A","-a"):
            pa = float(a);
            if pa > 1 or pa < 0:
                print "probability should between 0 to 1"
                sys.exit(2)
        elif o in ("-C","-c"):
            pc = float(a);
            if pc > 1 or pc < 0:
                print "probability should between 0 to 1"
                sys.exit(2)
        elif o in ("-G","-g"):
            pg = float(a);
            if pg > 1 or pg < 0:
                print "probability should between 0 to 1"
                sys.exit(2)
        elif o in ("-T","-t"):
            pt = float(a);
            if pt > 1 or pt < 0:
                print "probability should between 0 to 1"
                sys.exit(2)
        elif o == "-N":
            pn = float(a);
            if pn > 1 or pn < 0:
                print "probability should between 0 to 1"
                sys.exit(2)
        elif o in ("-h","--help"):
            print usage()
            sys.exit(0)
        elif o =="-v":
            coverage = int(a,10)
            if coverage < 0:
                print " the coverage of base paire should be non negative"
                sys.exit(0)
        elif o == "-k":
            lenOfkmer = int(a)
            if lenOfkmer < 0:
                print " the length of kmer should be non negative"
                sys.exit(0)
        elif o == "-l":
            lenOfreads = int(a)
            if lenOfreads < 0:
                print " the length of reads should be non negative"
                sys.exit(0)
        elif o == "-o":
            prefix = a
        elif o=="-s":
            lenOfseq = int(a,10)
            if lenOfseq < 0:
                print " the length of reads should be non negative"
                sys.exit(0)
        else:
            usage()
            sys.exit(0)
    # check basic conditions
    # pdb.set_trace()
    p = pa + pc + pg + pt + pn
    if p != 1:
        print 'The total probability of A({0}),C({1}),G({2}), T({3}) and N({4}) - {5} is not equal to one '.format(pa,pc,pg,pt,pn,p)
        sys.exit(0)
    if lenOfkmer > lenOfreads:
        raise RuntimeError( 'length of kmer({0}) should less or equal to the length of reads({1})'.format(lenOfkmer,lenOfreads))
        sys.exit(0)
    if lenOfseq == 0:
        print( 'length of whole sequence has to be specified')
        print "****************************usage******************************"
        usage()
        sys.exit(0)
    if lenOfseq < lenOfreads:
        print 'the length of whole sequence({0}) should be large \
than that of single read({1})'.format(lenOfseq,lenOfreads)
        sys.exit(0)
    print ' Arguments: probability of A({0}),\n \
           probability of C({1}), \n \
           probability of G({2}), \n \
           probability of T({3}), \n \
           probability of N({4}), \n \
           length of reads ({5}), \n \
           length of kmer({6}), \n \
           Coverage of base pair({7})\n \
           Length of whole sequence({8})\n \
           output prefix({9})'.\
           format(pa,\
                  pc,\
                  pg,\
                  pt,\
                  pn,\
                  lenOfreads,\
                  lenOfkmer,\
                  coverage,\
                  lenOfseq,\
                  prefix)
    run(pa,pc,pg,pt,pn,lenOfreads,lenOfkmer,coverage,lenOfseq,prefix)
if __name__=="__main__":
    driver()

