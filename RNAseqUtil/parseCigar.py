#!/usr/bin/pyhton
import re
def generateAlignmentSeq(cigar,seq):
    cut = re.split('M|I|D',cigar)[0:-1]
    off = 0
    roff = 0
    moff = 0
    editseq = ''
    for p in cut:
        off = off + len(p)
        if cigar[off] == 'M':
            editseq = editseq + seq[roff:roff+int(p)]
            roff = roff + int(p)
            moff = moff + int(p)
            off = off + 1
        elif cigar[off] == 'I':
            roff = roff + int(p)
            off = off + 1
        elif cigar[off] == 'D':
            editseq = editseq + int(p)*' '
            off = off + 1
            moff = moff + int(p)
    return editseq

