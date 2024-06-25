import itertools
import numpy as np

from enum import Enum
from ._maths import quat2dec

class Chm(Enum):
    chr1 = 1
    chr2 = 2
    chr3 = 3
    chr4 = 4
    chr5 = 5
    chr6 = 6
    chr7 = 7
    chr8 = 8
    chr9 = 9
    chr10 = 10
    chr11 = 11
    chr12 = 12
    chr13 = 13
    chr14 = 14
    chr15 = 15
    chr16 = 16
    chr17 = 17
    chr18 = 18
    chr19 = 19
    chr20 = 20
    chr21 = 21
    chr22 = 22
    chrX  = 23
    # chrY  = 24
    # chrM  = 25
    
BigWigChromSizesDict = {
    Chm(1): 249250621,
    Chm(2): 243199373,
    Chm(3): 198022430,
    Chm(4): 191154276,
    Chm(5): 180915260,
    Chm(6): 171115067,
    Chm(7): 159138663,
    Chm(8): 146364022,
    Chm(9): 141213431,
    Chm(10): 135534747,
    Chm(11): 135006516,
    Chm(12): 133851895,
    Chm(13): 115169878,
    Chm(14): 107349540,
    Chm(15): 102531392,
    Chm(16): 90354753,
    Chm(17): 81195210,
    Chm(18): 78077248,
    Chm(19): 59128983,
    Chm(20): 63025520,
    Chm(21): 48129895,
    Chm(22): 51304566,
    Chm(23): 155270560, # X
    # Chm(24): 59373566  # Y
    # Chm(25): 16571,    # M
}

class Nucl(Enum):
    A = 0
    C = 1
    G = 2
    T = 3

def build_N_gram_nucl(N):
    pwr = np.power(10, np.ones(N).cumsum()-1)[::-1] 
    elements = [ e.value for e in Nucl ] 
    ngram_dict = {}
    for ctx in itertools.product(elements, repeat=N):
        key = ''.join([ Nucl(c).name for c in ctx ])
        value = int(f"{sum(pwr*np.array(ctx))}")
        ngram_dict[key] = value

    return Enum("context3", ngram_dict, type=int)
