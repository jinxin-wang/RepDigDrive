import itertools
import numpy as np

from enum import Enum
# from ._maths import quat2dec

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

class NUCLEOTIDE_ENUM(Enum):
    C = 0
    T = 1
    A = 2
    G = 3

nucl = [ n.name for n in NUCLEOTIDE_ENUM ]

BASE_SUBSTITUTION_TYPES = {
    'C>A': 1, 'C>G': 2, 'C>T' : 3,
    'G>T': 1, 'G>C': 2, 'G>A' : 3,
    'T>A': 4, 'T>C': 5, 'T>G' : 6,
    'A>T': 4, 'A>G': 5, 'A>C' : 6,
    'INDEL': 0
}

BASE_SUBSTITUTION_CLASSES = {
    'C>A': 1, 'G>T': 1,
    'C>G': 2, 'G>C': 2, 
    'C>T': 3, 'G>A': 3, 
    'T>A': 1, 'A>T': 1, 
    'T>C': 3, 'A>G': 3, 
    'T>G': 2, 'A>C': 2,
    'INDEL': 0
}

class COHORT_CLASS(Enum):
    # TODO
    pass

class MUT_ANNOT(Enum):
    INDEL = 0
    Missense   = 1
    Nonsense   = 2
    Noncoding  = 3
    Stop_loss  = 4
    Synonymous = 5
    Essential_Splice = 6

def build_N_gram_nucl_enum(N):
    """
    for example : 3-gram 
    AAA: int("222", base = 4) 
    TCG: int("103", base = 4) 
    TTT: int("111", base = 4) 

    AAA + TTT = int("333", base = 4) = 63
    CTG + GAC = int("333", base = 4) = 63
    """
    name = f"context{N}"

    pwr = np.power(10, np.ones(N).cumsum()-1)[::-1] 
    elements = [ e.value for e in NUCLEOTIDE_ENUM ] 
    ngram_dict = {} 
    for ctx in itertools.product(elements, repeat=N):
        key = ''.join([ NUCLEOTIDE_ENUM(c).name for c in ctx ])
        # convert Quaternary to Decimal
        value = int(f"{int(sum(pwr*np.array(ctx)))}", base = 4)
        ngram_dict[key] = value

    return Enum(name, ngram_dict, type=int)

@DeprecationWarning
def build_N_ctx_mut_enum(name, N):
    """
    for example: 3-gram

    TCA -> A     : int('1022', base = 4)
    TCA -> INDEL : int('1020', base = 4)
    TCA -> G     : int('1023', base = 4)
    TCA -> T     : int('1021', base = 4)

    TTA -> A     : int('1122', base = 4)
    TTA -> C     : int('1120', base = 4)
    TTA -> G     : int('1123', base = 4)
    TTA -> INDEL : int('1121', base = 4)

    AAT -> INDEL : int('2212', base = 4)
    AAT -> C     : int('2210', base = 4)
    AAT -> G     : int('2213', base = 4)
    AAT -> T     : int('2211', base = 4)

    'AAT -> T' + 'TTA -> A' = int('2211', base = 4) + int('1122', base = 4) = int('3333', base = 4)
    """
    mut_enum = build_N_gram_nucl_enum(name, N)

    INDEL = [ m for m in mut_enum if m.name[int(N/2)-1] == m.name[-1]]
    SNV   = [ m for m in mut_enum if m.name[int(N/2)-1] != m.name[-1]]

    return mut_enum, INDEL, SNV

