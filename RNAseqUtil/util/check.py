from util.directory.folder import *
from util.command.subproc import *
import subprocess

def finished():
	proc = subprocess.Popen( ["-c", " ps -u luzhao | grep 'java' | wc -l" ],stdout=subprocess.PIPE,shell=True)
	out = int( proc.communicate()[0])

	if out > 0:
		return False
	else:
		return True
