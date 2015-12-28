import shlex,subprocess

def execute(command,out,err):
	command_line = shlex.split(command)
	p = subprocess.Popen(command_line,stdout = out,stdin = err)

