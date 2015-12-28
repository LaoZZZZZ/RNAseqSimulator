#!usr/bin/python


# find the maximum and minimum longest common string
class PrefixSpan(object):
	@staticmethod
	def findMaxMin(ref,pattern,start = 0,end = float('inf')):
		stat = {}
		if pattern and ref:
			pattern = pattern.upper()
			if end > len(ref):
				end = len(ref)
			ref = ref.upper()
			pos = list(range(start,end))
			off = list([0]*(end - start))
			res = zip(pos,off)
			for c in pattern:

				tmp = PrefixSpan.findcommon(ref,c,res)
				# the candidate become unique for the first time
				if len(tmp) == 1 and not stat:
					stat['min'] = tmp[0]
					stat['max'] = tmp[0]
				# keep updating the max unique common string after the candidate appear
				elif len(tmp)== 1 and stat:
					stat['max'] = tmp[0]
			                # stop after there is no candidate there
		                elif not tmp and stat:
					break
		                res = tmp
			# there are multiple position for this pattern
			if not stat and res:
				stat['min'] = res[0]
				stat['max'] = res[0]
		return stat
	@staticmethod
	def findMaxMinEnd(ref,pattern,start=0,end=float('inf')):
		stat = {}
		if pattern and ref:
			pattern = pattern.upper()
			if end > len(ref):
				end = len(ref)
			ref = ref.upper()
			pos = list(range(start,end))[::-1]
			off = list([0]*(end - start))
			res = zip(pos,off)
			for c in reversed(pattern):
				tmp = PrefixSpan.findCommonEnd(ref,c,res)
				# the candidate become unique for the first time
				if len(tmp) == 1 and not stat:
					stat['min'] = (tmp[0][0]-tmp[0][1]+1,tmp[0][1])
					stat['max'] = (tmp[0][0]-tmp[0][1]+1,tmp[0][1])
				# keep updating the max unique common string after the candidate appear
				elif len(tmp)== 1 and stat:
					stat['max'] = (tmp[0][0]-tmp[0][1]+1,tmp[0][1])
				# stop after there is no candidate there
				elif not tmp and stat:
					break
				res = tmp
			# there are multiple position for this pattern
			if not stat and res:
				stat['min'] = (res[0][0]-res[0][1]+1,res[0][1])
				stat['max'] = (res[0][0]-res[0][1]+1,res[0][1])
		return stat
	@staticmethod
	def findcommon(ref,c,pos):
		res = []
		length = len(ref)
		for o,off in pos:
			if o + off < length:
				if ref[o+off] == c:
					res.append((o,off+1))

		return res
	@staticmethod
	def findCommonEnd(ref,c,pos):
		res = []
		for o,off in pos:
			if o - off >= 0:
				if ref[o - off] == c:
					res.append((o,off+ 1))
		return res
if __name__ == '__main__':

	ref = 'caaaggcagcaagaag'
	pat = 'caaaggaag'
	res = PrefixSpan.findMaxMin(ref,pat)
	if res:
		print(res['min'],res['max'])
		print(ref[res['min'][0]:res['min'][0]+res['min'][1]],ref[res['max'][0]:res['max'][0]+res['max'][1]])
