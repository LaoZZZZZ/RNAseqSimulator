#!/usr/bin/python

import numpy as np
from fragmentGenerator import fragmentGenerator as fGenerator

from mrnaDetermineGenerator import mrnaDetermineGenerator
from amplification import amplifier
from reads import read
from mrnalibrary import mrnaLibrary
from sequenceError import sequenceError
from sequencer import Sequencer
from sequenceError import sequenceError
from readsGenerator import readGenerator
from pairSequencer import pairSequencer
import getopt
from expression_stats import expression_stats
def rnaSimulator(template,rdlen,size,lambd,filter_lower,filter_upper,pair,outpref,seed):
	np.random.seed(int(seed))
	mrnalib = mrnaLibrary(template) #object that holds all template 

	mrnaGen = mrnaDetermineGenerator(template,int(size))	
	fragGen = fGenerator(mrnaGen,int(lambd),int(filter_lower),int(filter_upper))  # generate all fragments given the transcriptome
	
	ErrorModel = sequenceError()  # define the  sequence error model
	seqmachine = Sequencer(ErrorModel,int(rdlen),int(pair))
	
	simulator = readGenerator(seqmachine,mrnalib,fragGen,outpref)
	simulator.simulatorDrive() 
	expression_stats(mrnaGen.getExpressionlevels(),outpref+'.sam',outpref+'_expressionLevel.csv').dump()
def drive():
	setting = {}
	setting['template'] ='/home/luzhao/simulation/simulate_data/ref/mrna.fa'  #template library 
	setting['rdlen'] = 50
	
	setting['lambd'] = 120  # fragment average length
	setting['filter_lower'] = 100 # lower bound of fragment's length  
	setting['filter_upper'] = 140 # upper bound of fragment's length 
	setting['pair'] = True
	setting['seed'] = 23
	setting['size'] = 10
	setting['outpref'] = '/home/luzhao/simulation/simulate_data/simulated'
	rnaSimulator(**setting)

if __name__=='__main__':
	
	drive()
