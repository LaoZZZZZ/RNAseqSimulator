#!/usr/bin/python

import os
import sys
from extract_mapped import source_parsing
if __name__ == '__main__':
	task = source_parsing(sys.stdin)
	dest = os.path.join(task['dir'],task['targets'][0])
	fhand = open(dest,'r')
	total = 0
	for line in fhand:
		total = total + 1
		if total > 16:
			break
		print(line.split()[0])
