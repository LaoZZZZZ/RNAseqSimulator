#!/usr/bin/python

def simple_generator():
	for i in range(10):
		yield i
if __name__ == '__main__':
	for i in simple_generator():
		print(i)
