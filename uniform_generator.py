#!/usr/bin/python

#from numpy import *
import numpy as np


class uniformGenerator(object):
    def __init__(self,seed = 2147483647):
        np.random.seed(seed)
    def __iter__(self):
        return self
    def __next__(self):
        return np.random.uniform(0,1,1)



if __name__ == '__main__':
    unifGen = uniformGenerator()
    for i in unifGen:
        print(i)
        if i < 0.1:
            break
