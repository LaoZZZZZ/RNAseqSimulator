import os
import sys
# show all files in this directory and its subdirectory
def allFiles( directory):
	lst = []
	for root,dirs,files in os.walk(directory):
		for f in files:
			lst.append(f)
	return lst
#show all subdirectory in this dir

def allDirs(directory):
	lst=[]
	for root,dirs,files in os.walk(directory):
		for d in dirs:
			lst.append(d);
	return lst
# only show the current level of subdirectory
def currentDirs(directory):
	lst = []
	for root,dirs,files in os.walk(directory):
		for d in dirs:
			pd = os.path.join(directory,d)
			if os.path.isdir(pd):
				lst.append(d)
	return lst
# only show those files in this directory, no subdirectory is shown
def currentFiles(directory):
	lst = []
	for root,dirs,files in os.walk(directory):
		for f in files:
			pf = os.path.join(directory,f)
			if os.path.isfile(pf):
				lst.append(f)
	return lst

#show all files and subdirectory in this dir

def ls(directory):
	return currentFiles(directory) + currentDirs(directory)
#check the existence of  the file in the directory

def existfile(dir,file):
	lst = currentFiles(dir)
	return file in lst
def existfilethorough(dir,file):
	lst = allFiles(dir)
	return file in lst
# check the existence of the specific direcoty
def existdir(dir,subdir):
	return subdir in currentDirs(dir)
def existdirthorough(dir,subdir):
	return subdir in allDirs(dir) 

