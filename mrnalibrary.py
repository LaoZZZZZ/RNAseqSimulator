#!/usr/bin/python


# mrna library that contains all transcript sequence
class mrnaLibrary(object):
	def __init__(self,mrnafile):
		self.loadMrna(mrnafile)
	def loadMrna(self,mrnafile):
		data = open(mrnafile,'r').read().split('\n')
		i = 0 
		self.mrnas = []
		self.mrnaDict = {}
		while(i < len(data)-1 ):
			self.mrnas.append((data[i].strip('>'),data[i+1]))
			self.mrnaDict[data[i].strip('>')] = i/2 
			i += 2
	def __iter__(self):
		return self
	def __getitem__(self,idx):
		if isinstance(idx,basestring):
			return self.mrnas[self.mrnaDict[idx]] 
		elif idx >= len(self.mrnas):
			raise IndexError('out of range')
		return self.mrnas[idx]

	def __len__(self):
		return len(self.mrnas)
	#def getMrnas(self):
	#	return dict(self.mrnas)



if __name__ == '__main__':
    mrnafile = 'mrna.txt'
    lib = mrnaLibrary(mrnafile)
    print(lib['why'],lib[1])
        
            
        
