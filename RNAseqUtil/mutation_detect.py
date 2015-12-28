#!/usr/bin/python
import os
import sys
from alignment import alignment
from align_candidate import align_candidate
from prefixspan import PrefixSpan
from reads import read
from mutation import mutation
from strMatch import *
from alignment import alignment
from reads import read
def trim(id):
        tmp = id.split('/')[0]
        if tmp[0] == '@':
                tmp = tmp[1:]
        return tmp

def loadSamImp(samfile,start = 0,end = float('Inf')):
	print(samfile)
	dictsam = {}
	if not os.path.isfile(samfile):
		print('{0} is not a sam file'.format(samfile))
		sys.exit()
	fhandle = open(samfile,'r')
	for line in fhandle:
		if line[0] != '@':
			rec = line.rstrip('\n').split()
			if len(rec) > 9:
				startpos = int(rec[3])
				if startpos >= start and startpos <= end:
					al = alignment(trim(rec[0]),rec[2],int(rec[1]),startpos,rec[9],rec[5])
					if al.getID() in dictsam:
						dictsam[al.getID()].append(al)
					else:
						dictsam[al.getID()] = [al]
	return dictsam
def loadPartialSamByPos(samfile,mrnas):
	print(samfile)
	dictsam = {}
	if not os.path.isfile(samfile):
		print('can not find the file %s' %(samfile))
	fhandler = open(samfile,'r')
	for line in fhandler:
		if line[0] != '@':
			rec = line.rstrip('\n').split()
			if not rec[2] in mrnas:
				continue
			al = alignment(rec[0],rec[2],rec[1],rec[3],rec[9],rec[5])
			if rec[2] in dictsam:
				dictsam[rec[2]].append(al)
			else:
				dictsam[rec[2]] = [al]
	return dictsam
def loadSamByPos(samfile):
	print(samfile)
	dictsam = {}
	if not os.path.isfile(samfile):
		print('can not find the file %s' %(samfile))
	fhandler = open(samfile,'r')
	for line in fhandler:
		if line[0] != '@':
			rec = line.rstrip('\n').split()
			al = alignment(rec[0],rec[2],rec[1],rec[3],rec[9],rec[5])
			if rec[2] in dictsam:
				dictsam[rec[2]].append(al)
			else:
				dictsam[rec[2]] = [al]
	return dictsam
def loadPairSamImp(samfile):
	print(samfile)
	dictsam = {}
	if not os.path.isfile(samfile):
		print('{0} is not a sam file'.format(samfile))
		sys.exit()
	fhandle = open(samfile,'r')
	for line in fhandle:
		if line[0] != '@':
			rec = line.rstrip('\n').split()
			if len(rec) > 9:
				if rec[5] != '*':
					startpos = int(rec[3])
					al = alignment(trim(rec[0]),rec[2],int(rec[1]) & 0x10,startpos,rec[9],rec[5])
					if rec[2] in dictsam:
						if al.getID() in dictsam[rec[2]]:
							dictsam[rec[2]][al.getID()].append(al)
						else:
							dictsam[rec[2]][al.getID()] = [al]
					else:
						dictsam[rec[2]] = {al.getID():[al]}
	return dictsam
def loadPairSam(samdir,target):
	import os
	samfile = os.path.join(samdir,target+'.sam')
	return loadPairSamImp(samfile)
def loadSam(samdir,target,start = 0,end = float('Inf')):
	import os
        samfile = os.path.join(samdir,target+'.sam')
       	return loadSamImp(samfile,start,end) 
def loadfqimpl(fqfile):
        #print(fqfile)
        dictfq = {} 
	import os
	from reads import read
        if not os.path.isfile(fqfile):
                print('{0} is not a fq file'.format(fqfile))
                sys.exit()
        fhandle = open(fqfile,'r')
        rd = []
        total = 0
        for line in fhandle:
		total = total + 1
		rd.append(line.rstrip('\n'))
		if total % 4 == 0:
			tmp = read(trim(rd[0]),rd[1].upper(),rd[3])
			if tmp.getid() in dictfq:
				print('there are duplicate reads %s \n'%(rd[0]))
			else:
                        	dictfq[tmp.getid()] = tmp
                        rd = []
        return dictfq

def loadfq(fqdir,target):
        fqfile = os.path.join(fqdir,target)
        return loadfqimpl(fqfile)
def filter(sam,fqdir,target):
    fq = os.path.join(fqdir,target+'.fq')
    fhandler = open(fq,'r')
    output =  open(os.path.join(fqdir,target+'region.fq'),'w')
    rd = []
    total = 0
    for line in fhandler:
        total = total + 1
        rd.append(line.rstrip('\n'))
        if total % 4 == 0:
            if trim(rd[0]) in sam:
                output.write('%s\n%s\n+\n%s\n'%(rd[0],rd[1],rd[3]))
            rd = []
    output.close()
def revcompl(seq,dictm = {'A':'T','T':'A','C':'G','G':'C','N':'N'}):
        tm = ''
        for e in seq:
            tm = dictm[e] + tm
        return tm
def findCommon(seq,alignments):
        head = []
        tail = []
        if alignments:
                length = len(alignments[0].getSeq()) 
                for e in alignments:
                        if e.getOrient() == 0:

                                if findMismatches(seq[0:length],e.getSeq()) <= 1:
                                        head.append(e)
                                else:
                                        #print(seq[len(seq)-length:],e)
                                        if findExactMismatches(seq[len(seq)-length:],e.getSeq()) <= 1:
                                        	tail.append(e)
					else:
						print(seq[len(seq)-length:],e.__str__())
                        else:
                                
                                if findMismatches(seq[0:length],revcompl(e.getSeq())) <= 1:
                                        head.append(e)
                                else:
                                        if findExactMismatches(seq[len(seq)-length:],revcompl(e.getSeq())) <=1: 
                                        	tail.append(e)
					else:
						print(seq[len(seq)-length:],e.__str__())
        return [head,tail]

#find those valid match between head and tail alignment based on the original reads and the position information
def findmatch(heads,tails,rd):
        match = []
        if heads and tails:
                seqlength = len(heads[0])
                for h in heads:
                        for t in tails:
                                if h.getOrient() == t.getOrient():
                                        if h.getPos() + seqlength <= t.getPos():
                                                match.append((h,t))

        elif heads and not tails:
                for h in heads:
                        match.append((h,'na'))
        elif not heads and tails:
                for t in tails:
                        match.append(('na',t))
        return match
# validate all the alignements of one read based on the alignment information and original read sequence                        
def validateSingle(rd,alignments):
        cand = align_candidate(rd)
        if alignments:
                [heads,tails] = findCommon(rd.getseq().upper(),alignments)
                match = findmatch(heads,tails,rd)
                cand.addmatch(match)
        return cand
        
# validate all the alignement based on the alignment information and original reads sequence
def validate(alignments,readsRec):
        matches = []
        if alignments and readsRec:
                for id,al in alignments.items():
                        rd = readsRec[id]
                        if rd:
                                res = validateSingle(rd,al)
                                if len(res):
                                        matches.append(res)
        return matches
#def loadRef(refdir,target):

                        

# extent the alignment along the direction specified, find the maximum valid alignment 
def extendAlignment(align,rdseq,ref,head = True):
        from alignment import alignment
        match = ''
        pos = align.getPos()
        ref = ref.upper()
        rdseq = rdseq.upper()
        #print(rdseq)
        if head:
                pos = pos + len(align)
                l = len(ref)
                left = rdseq[len(align):]
                for p in range(len(left)):
                        if l > pos + p - 1 and left[p] == ref[pos + p - 1]:
                                match = match + left[p]
                        else:
                                break
                #print(match)
                #print(align)
                return alignment(align.getID(),align.getChr(),align.getOrient(),align.getPos(),align.getSeq()+match)
                #print(r,align)
                #return r
        else:
    
                left = rdseq[0:len(rdseq) - len(align)]
                l = len(left)
                #print(left,len(left),pos,pos-len(left),ref[pos-len(left)-1:pos-1])
                #print(align)
                for p in range(l):
                        #print(left[-p-1],ref[pos - p - 2],ref[pos-1])
			#if pos - p -2 >= len(ref):
				#print(len(ref),pos-p-2,p,)
                        if 0 <= pos - p - 2 and pos -p -2 < len(ref) and left[-p-1] == ref[pos - p - 2]:
                                match = left[-p - 1] + match
                        else:
                                break
                #print(match)
                
                return alignment(align.getID(),align.getChr(),align.getOrient(),pos - len(match),match+align.getSeq(),str(len(match+align.getSeq())) +'M') 
                #print(r,align)
                #return r
                        
      
# 0 is normal, 1 is insertion, 2 is deletion,3 is snp
def specifyMutations(head,tail,rd,ref,rdseq):
        try:
                from mutation import mutation
                from align_candidate import align_candidate
                lhead = len(head)
                ltail = len(tail)
                lrd = len(rd)
                mutationPos = 0
                if lhead + ltail == lrd:
                        if head.getPos() + lhead == tail.getPos():
                                return None
                        if head.getPos() + lhead < tail.getPos():
                                #print(head,tail)
                                m = align_candidate(rd)
                                m.addmatch([(head,tail)])
                                return mutation(2,head.getChr(),head.getPos() + lhead,ref[head.getPos() + lhead-1:tail.getPos()-1],[m])
                        # there is overlap between head and tail
                        else:
                                overlaph = head.getSeq()[tail.getPos()-head.getPos():lhead]
                                overlapt = tail.getSeq()[0:head.getPos() + lhead - tail.getPos()]
                                if overlaph == overlapt:
                                        m = align_candidate(rd)
                                        m.addmatch([(head,tail)])
                                        return mutation(1,head.getChr(),head.getPos() + lhead,overlaph,[m])
                                        #return mut
                                else:
                                        return None
                
                elif lhead + ltail < lrd:
                        if head.getPos() + lhead == tail.getPos():
                                #print(lhead,ltail,len(rd),rdseq[lhead:lrd - ltail],head,tail,head.getPos() + lhead)
                                m = align_candidate(rd)
                                m.addmatch([(head,tail)])
                                return mutation(1,head.getChr(),head.getPos() + lhead,rdseq[lhead:lrd - ltail],[m])
                        else:
                                return None
                else:
                        raise Exception('invalid head(%s) and tail(%s) alignment'%(head.__str__(),tail.__str__()))
        except Exception as err:
                print(err)
                import sys
                sys.exit(0)
                
                                        
def validateMutations(ref,matches):
        from prefixspan import PrefixSpan
        from mutation import mutation
        from alignment import alignment
        start = 0
        end = float('inf')
        result = []
        #print(matches)
        for head,tail in matches.getMatch():
                if head != 'na':
                                     
                        # both head and tail could be mapped to the reference
                        if tail != 'na':
                        
                                assert(head.getOrient()==tail.getOrient())
                                assert(head.getPos() < tail.getPos())
                                if head.getOrient() == 16:
                                        rdseq =  revcompl(matches.getRead().getseq())
                                else:
                                        rdseq = matches.getRead().getseq()
                                pat = rdseq[len(head):len(matches.getRead())-len(tail)]
                                
                                #pat = rdseq[0:len(matches.getRead())-len(tail)]
                                start = head.getPos() + len(head)
                                end = tail.getPos()
                                assert(start <= end)         
                                # there exists a possible insertion
                                
                                if end - start != len(rdseq) - len(head) - len(tail):
                                        # if the head and tail are adjacent, then this is an insertion
                                        rd = head.getSeq() + pat
                                        tmphead = extendAlignment(head,rd,ref)
                                        mut = specifyMutations(tmphead,tail,matches.getRead(),ref,rdseq)
                                        if not mut:
                                                tmptail = extendAlignment(tail,rd[len(tmphead):] + tail.getSeq(),ref,False)
                                                mut = specifyMutations(tmphead,tmptail,matches.getRead(),ref,rdseq)
                                                if mut:
                                                        result.append(mut) 
                                        else:
                                                result.append(mut)
                                else :
                                        # if the mismatch between the middle part of read and the ref smaller than 1, then report it
                                        # report possible mutation
                                        #print(ref[head.getPos() + len(head)-1:tail.getPos()-1].upper(),pat.upper())
                                        if findExactMismatches(ref[head.getPos() + len(head)-1:tail.getPos()-1].upper(),pat.upper()) == 1:
                                                allel = ''
                                                pos = 0
                                                for p in range(len(pat)):
                                                        if head.getPos() + len(head)+p-1 < len(ref) and pat[p].upper() != ref[head.getPos() + len(head)+p-1].upper():
                                                                allel = pat[p]
                                                                pos = head.getPos() + len(head) + p - 1
                                                                break
                                                head = alignment(head.getID(),head.getChr(),head.getOrient(),head.getPos(),head.getSeq()+pat)
                                                m = align_candidate(matches.getRead())
                                                m.addmatch([(head,tail)])
                                                mut = mutation(3,head.getChr(),pos,allel,[m])
                                                result.append(mut)
                        else:
                                         
                                 if head.getOrient() == 16:
                                         rdseq =  revcompl(matches.getRead().getseq())
                                 else:
                                         rdseq = matches.getRead().getseq()
                                 # first extend the head alignment as much as possible
                                 tmphead = extendAlignment(head,rdseq,ref)
                                 if len(tmphead) < len(matches.getRead()):
                                         
                                 #result.append(head)
                                         # extract the unmapped pat
                                         pat = rdseq[len(tmphead):]
                                         # search from tail to find the maximum match position
                                         tailpos = PrefixSpan.findMaxMinEnd(ref,pat)
                                         # only consider those tail part larger than 10 bp
                                         if tailpos and tailpos['max'][1] > 10:
                                                 #print(head,tmphead,len(tmphead),len(matches.getRead()))

                                                 tailpos = tailpos['max']
                                                 tmptail = alignment(head.getID(),head.getChr(),head.getOrient(),tailpos[0] + 1,ref[tailpos[0]:tailpos[0]+tailpos[1]])
                                                 #print(tailpos,tmptail)

                                                 mut = specifyMutations(tmphead,tmptail,matches.getRead(),ref,rdseq)
                                                 #print(mut)
                                                 if mut:
                                                         result.append(mut)
                                         else:
                                                 # if the last four base paire are different, then report this event
                                                 if len(pat) < 4:
                                                         m = align_candidate(matches.getRead())
                                                         m.addmatch([(tmphead,None)])
                                                         mut = mutation(1,tmphead.getChr(),len(tmphead) + tmphead.getPos(),pat,[m])
                                                         result.append(mut)       

                else:
                        if tail == 'na':
                                pass
                        else:
                                if tail.getOrient() == 16:
                                        rdseq =  revcompl(matches.getRead().getseq())
                                else:
                                        rdseq = matches.getRead().getseq()
                                tmptail = extendAlignment(tail,rdseq,ref,False)
                                if len(tmptail) < len(matches.getRead()):
                                        pat = rdseq[0:len(rdseq) - len(tmptail)]
                                        headpos = PrefixSpan.findMaxMin(ref,pat)
                                        if headpos and headpos['max'][1] > 10:
                                                headpos = headpos['max']
                                                tmphead = alignment(tail.getID(),tail.getChr(),tail.getOrient(),headpos[0]+1,ref[headpos[0]:headpos[0]+headpos[1]])
                                                mut = specifyMutations(tmphead,tmptail,matches.getRead(),ref,rdseq)
                                                if mut:
                                                        result.append(mut)
                                        else:
                                                # if the last four base aire are different, then report this event
                                                if len(pat) < 4:
                                                        m = align_candidate(matches.getRead())
                                                        m.addmatch([(None,tmptail)])
                                                        mut = mutation(1,tmptail.getChr(),tmptail.getPos() - len(pat),pat,[m])
                                                        result.append(mut)   
        if len(result) > 1:
                print('the read (%s) supports multiple mutation sigunature, which should be examined futher'%(matches.getRead().__str__()))
        return result
def drive(dictfq,samfile,reffile,start,end,outdir,supportCutOff = 2):
	from refKmer import readfa
        aligns = loadSamImp(samfile,start,end)
        #print(len(aligns))
        #filter(aligns,inputdir,target)
       	ref = readfa(reffile) 
        result = validate(aligns,dictfq)
        #print(len(result))
        mutations = {}
        for match in result:
                #final.append(realign1(ref['seq'][0],match))
                tmp = validateMutations(ref['seq'][0],match)
                #print(tmp[0])
                for m in tmp:
                        if m.__hash__() in mutations:
                                mutations[m.__hash__()].addSupport(m.getSupport())
                        else:
                                mutations[m.__hash__()] = m
      	support = ''
	mutfile = ''
	for k,v in mutations.items():
		su = v.getSupport()
		if len(su) >= supportCutOff:
			mutfile = mutfile + v.__str__() + '\n'
			for e in su:
				support = support + e.__rec__() + '\n'
	if support: 
		basename = os.path.splitext(os.path.basename(samfile))[0] 
		supp = open(os.path.join(outdir,basename+'.supp'),'w')
		supp.write(support)
		supp.close()
		mut = open(os.path.join(outdir,basename+'.mut'),'w')
		mut.write(mutfile)
		mut.close()
def getMappedCoverage(mut,dictsam,rdlen = 50):
	#start,end = mut.getChrMutRange()
	mrna=mut.mutationChr()
	mutpos = mut.mutationPos()
	mutlen = len(mut) 
	support = mut.getSupport()
	suppaligns = []
	for e in support:
		match = e.getMatch()
		for h,t in match:
			if isinstance(h,alignment):
				suppaligns.append(h)
			if isinstance(t,alignment):
				suppaligns.append(t)
	start = end = mutpos
	if mut.mutationType() == 1 or mut.mutationType() == 3:
		start = mutpos - rdlen + 1
	elif mut.mutationType() == 2:
		start = mutpos - rdlen + 1
		end = mutpos +  mutlen - 1
	over_aligns = []
	aligns = []
	if mrna in dictsam:
		aligns = dictsam[mrna]
	for a in aligns:
		if a.getPos() >= start and a.getPos() <= end:
			over_aligns.append(a) 
	#while start <= end:
	#	if mrna +':' + str(start) in dictsam:
	#		aligns = dictsam[mrna+':'+str(start)]
	#		over_aligns = over_aligns + aligns 
	#	start = start + 1
	res = []
	for al in over_aligns:
		if not al in suppaligns:
			rd = read(al.getID(),al.getSeq(),len(al.getSeq())*'#')
			match = align_candidate(rd)
			match.addmatch([(al,'')])
			res.append(match)
	return res	
def pairdrive(dictfq,splitsamfile,orig_samfile,ref,refgenes,outdir,supportCutOff = 2):
	from parseRefGene import loadRefMrna 
        aligns = loadPairSamImp(splitsamfile)
        #filter(aligns,inputdir,target)
       	#ref = loadRefMrna(reffile)
	result = {}
	# organize the alignment by the mrna id
	for k,v in aligns.items(): 
        	result[k] = validate(v,dictfq)
        #print(len(result))
        mutations = {}
	# seach the mutations for each rna
        for k,rmut in result.items():
		for match in rmut:
                	tmp = validateMutations(ref[k][0],match)
                	for m in tmp:
                       		if m.__hash__() in mutations:
                               		mutations[m.__hash__()].addSupport(m.getSupport())
                        	else:
                               		mutations[m.__hash__()] = m
	validmutation = []
      	support = ''
	mutfile = ''
	if mutations:
		dictsam = loadSamByPos(orig_samfile)
		#print(len(dictsam))
	# filter out those valid mutations that supported by at least supportCufOff unique reads
	for k,m in mutations.items():
	#	print(m.NumOfSupportReads())
		if m.NumOfSupportReads() >= supportCutOff:
			if not m.mutationChr() in ref:
                                raise('can not find mrna id %s in refgene file'%(m.mutationChr()))
			for e in refgenes[m.mutationChr().split('ex')[0]]:
				if e.rnaname() == m.mutationChr():
					if e.inCodingRegion(m.mutationPos()):
						m.setChromosome(e.chr())
						m.setChrMutRange(e)
						#m.setChrPos(e.getChrPos(m.mutationStartPos()))
						m.setGene(e.genename())
						normalsupport = getMappedCoverage(m,dictsam,50)
						m.addNormalSupport(normalsupport)
					#	m.setUniqMapp(uniq)
					#	m.setOverallMapp(total)	
						validmutation.append(m)
	validmutation.sort(reverse=True)
	#print(len(validmutation))
	for m in validmutation:
		mutfile = mutfile + m.__str__() + '\n'
		for e in m.getSupport():
			support = support + e.__rec__() + '\n'
	if mutfile and support:
		 
		basename = os.path.splitext(os.path.basename(orig_samfile))[0] 
		supp = open(os.path.join(outdir,basename+'.supp'),'w')
		supp.write(support)
		supp.close()
		mut = open(os.path.join(outdir,basename+'.csv'),'w')
		head = '@ 1:insertion \t 2:deletion \t 3:snp\n'
		head =  head + '@column name\n'
		head = head + '@mutation_type,mrna,mrnapos,chr,chr_mutation_start,chr_mtation_end,gene,mutation_length,mutation_seq,uniquesupport,totalsupport,uniquenorm,totalnorm,uniquepercent,totalpercent,supportReadsID,normalReadsID\n'
		mut.write(head)
		mut.write(mutfile)
		mut.close()
	 
if __name__ == '__main__':
	import os
    	import sys
	from refKmer import readfa
    	from alignment import alignment
    	from reads import read
    	from align_candidate import align_candidate
	from samDrive import splitMrna
	from parseRefGene import parseRefGene
	from refGeneStruct import refGene
    	directory = '/home/luzhao/THR104_analysis'
	refgenes = parseRefGene(os.path.join(directory,'chromo/refGene.txt'))
	mrnas = []
	for k,v in refgenes.items():
		for m in v:
			mrnas.append(m.rnaname())
	splits = splitMrna(mrnas,2)

    	fqfile = os.path.join(directory,'unmapped/THR104_unmapped.fq')
    	ref = os.path.join(directory,'chromo/ref/mrna.fa')
	orig_samfile = os.path.join(directory,'samfile/test.sam')
	for e in splits:
		dictsam = loadPartialSamByPos('/home/luzhao/THR104_analysis/mapped/THR104.sam',set(e))
		dictsam = {}
		#print(len(e))
#	total = 0
#	for k,v in dictsam.items():
#		total = total + v
#	print(total,len(dictsam))	
   	#print(ref['seq'][0][1168:1180],ref['seq'][0][1138:1160])
#    	start = 0
#    	end = float('inf')
#	outdir = directory 
#       dictfq = loadfqimpl(fqfile)
#	aligns = loadPairSam(directory,'THR104_split')
#	print(len(aligns))
#	from parseRefGene import loadRefMrna 
#	ref = loadRefMrna(ref)
#	splitsamfile = '/home/luzhao/THR104_analysis/test.sam'
#	from parseRefGene import parseRefGene
#	pairdrive(dictfq,splitsamfile,orig_samfile,ref,refgenes,outdir,1)
#	drive(dictfq,samfile,ref,start,end,outdir)

