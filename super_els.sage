from sage.all import *

def superelliptic_els(coefflist, k=2):
    raise NotImplementedError("This isn't finished yet!")
    if not R_soluble(coefflist, k):
        return False
    for p in ZZ(k).prime_divisors():
        if not Qp_soluble(coefflist, k):
            return False
    R.<x> = QQ[]
    