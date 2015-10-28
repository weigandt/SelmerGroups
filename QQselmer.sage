import sys
from sage.all import *

def QQ_selmer(S = [], m = 2):
#This is the analog of Robert Miller' K.selmer_group() for K == QQ.
    gens = [-1]
    for p in S:
        if not is_prime(p):
            raise RuntimeError("Element " + str(p) + " is not a prime!")
        if p in gens:
            raise RuntimeError("List S must be of DISTINCT primes!")
        else:
            gens.append(p)
    return gens