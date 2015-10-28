import sys
from sage.all import *

def four_selmer_rank_bound(E, remove_torsion = True):
    coeff = str(list(E.a_invariants()))
    magma.eval('E := EllipticCurve(' + coeff + ')')
    addin = ''
    if remove_torsion:
        addin = ': RemoveTorsion'
    magma.eval('Sel2,Sel2_map := TwoSelmerGroup(E ' + addin + ')')
    magma.eval('two_descent_basis := [TwoCover(Sel2.i @@ Sel2_map) : i in [1..Ngens(Sel2)]]')
    magma.eval('CTmatrix := Matrix(GF(2), Ngens(Sel2), [CasselsTatePairing(C,D) : C, D in two_descent_basis])')
    magma.eval('assert IsSymmetric(CTmatrix)')
    return ZZ(magma.eval('Ngens(Sel2)')) - ZZ(magma.eval('Rank(CTmatrix)'))

def p_torsion_rank(E,p):
    r"""
    Return the dimension of the p-torsion subgroup of 'E(K)'.
    
    This will be 0, 1  or 2.
    """
    if not p.is_prime():
        raise NotImplementedError("n_torsion_rank is not implemented for composite numbers.")
    elif p == 2:
        return E.two_torsion_rank()
    else:
        f = E.division_polynomial(p)
        n = 2*len(f.roots())+1
        return ZZ(n).ord(p)

def four_selmer_kernel_size(E):
    if E.two_torsion_rank() == 0:
        return ZZ(1)
    elif E.two_torsion_rank() == 1:
        if E.torsion_order() % 4 == 0:
            return ZZ(1)
        else:
            return ZZ(2)
    else:
        structure = E.torsion_subgroup().invariants() #This won't work for cyclotomic fields, etc. Use div. polynomials.
        four_torsion_rank = 0
        for n in structure:
            if n % 4 == 0:
                four_torsion_rank+=1
        return ZZ(2^(2 - four_torsion_rank))

def four_selmer_size(E):
    coeff = str(list(E.a_invariants()))
    magma.eval('E := EllipticCurve(' + coeff + ')')
    magma.eval('Sel2,Sel2_map := TwoSelmerGroup(E)')
    magma.eval('two_descent_basis := [TwoCover(Sel2.i @@ Sel2_map) : i in [1..Ngens(Sel2)]]')
    magma.eval('CTmatrix := Matrix(GF(2), Ngens(Sel2), [CasselsTatePairing(C,D) : C, D in two_descent_basis])')
    scoker = ZZ(magma.eval('Rank(CTmatrix)'))
    #print scoker
    s2 = ZZ(magma.eval('Ngens(Sel2)'))
    #print s2
    sker = four_selmer_kernel_size(E).ord(2)
    #print sker
    sel4exp = 2*s2 - scoker - sker
    #print sel4exp
    return 2^sel4exp
