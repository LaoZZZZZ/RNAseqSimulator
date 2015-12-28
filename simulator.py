import numpy as np
from reads import rdGenerator as rd
from seqGenerator import seqGenerator as seqgen 
def run(proba,probc,probg,probt,probn,lenOfreads,lenOfkmer,coverage,lenOfseq,prefix):
    generator = seqgen(proba,probc,probg,probt,probn,prefix,lenOfseq)
    # generate the sequence that has the specified length
    seq = generator.generateStraight()
    reads = rd(seq,coverage,lenOfreads,lenOfkmer,prefix)
    reads.generate()
    #print seq
    
