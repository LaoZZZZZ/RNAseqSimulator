#!/usr/bin/python
import numpy as np
from discrete_uniform_generator import discreteUniform
import math
class fragmentLength:

    # return the start position in the transcript
    @staticmethod
    def startPos(start,end):
        return discreteUniform.next(start,end,1)[0]
    # return the length of the fragment given the lambda
    # lambd is the expected length of the fragment, which is 1/lambda in the expo
    @staticmethod
    def endPos(lambd):
        obs = np.random.exponential(lambd,1)
        if obs - int(obs) > 0.5:
            return math.ceil(obs)
        else:
            return math.floor(obs)
    # return the start and length of the fragment
    @staticmethod
    def getFragment(start,end,lambd):
        return (fragmentLength.startPos(start,end),fragmentLength.endPos(lambd))



if __name__ == '__main__':
    print(fragmentLength.getFragment(0,1000,500))
