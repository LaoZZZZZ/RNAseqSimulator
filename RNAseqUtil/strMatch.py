#!/usr/bin/python
from prefixspan import PrefixSpan
def countMatch(ref,pat):
        from prefixspan import PrefixSpan
        result = PrefixSpan.findMaxMin(ref,pat)
        matched = 0
        while result:
                left1 = ref
                left2 = pat
                while result:
                        matched = matched + result['max'][1]
                        left1 = removeParts(left1,result['max'][0],result['max'][1])
                        left2 = left2[result['max'][1]:]
                        result = PrefixSpan.findMaxMin(left1,left2)
        return matched
def countExactMatch(ref,pat):
        try:
                if len(ref) < len(pat):
                        raise Exception('invalid comparison between the ref and pattern where the length of ref is smaller than that of patter')
                m = 0
                for pos in range(len(pat)):
                        if ref[pos] == pat[pos] or ref[pos] == 'N' or pat[pos] == 'N':
                                m = m + 1
                return m
        except Exception as err:
                print(err)
                import sys
                sys.exit(0)
def findExactMismatches(seq1,seq2):
        if len(seq1) > len(seq2):
                return len(seq1) - countExactMatch(seq1[0:len(seq2)],seq2)
        else:
                return len(seq2) - countExactMatch(seq2[0:len(seq1)],seq1)
def findMismatches(seq1,seq2):
        from prefixspan import PrefixSpan
        matched = 0
        if len(seq1) > len(seq2):
                return len(seq1) - countMatch(seq1,seq2)
        else:
                return len(seq2) - countMatch(seq2,seq1)

def removeParts(seq,start,off):
        try:
                if len(seq) < off or len(seq) < start + off:
                        raise Exception( 'invalid range to remove ')
                if len(seq) > start + off:
                        back = seq[start+off:]
                        return back
                else:
                        return ''
        except Exception as err:
                print(err)
                import sys
                sys.exit(0)



