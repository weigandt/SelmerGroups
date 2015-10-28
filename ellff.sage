import sys
from sage.all import *

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

def quartic_to_jacobian(T):
    r"""
    INPUT: Coefficients (a0,a1,a2,a3,a4) of a binary quartic form. 
    OUTPUT: The Jacobian of the quartic curve C: y^2 = a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0 as an elliptic curve E, together with a rational map from C to E.
    """
    raise NotImplementedError('You need to get to work on this!')

def frey_to_sw(A,C):
    r"""
    INPUT: Coeffiecents A,C the Frey curve E: y^2 = x(x+A)(x - C)
    OUTPUT: Coefficents A_short and B_short of a short weierstrass equation for E together with an isomorphism mapping the Frey model to the short weiertrass model. 
    """
    raise NotImplementedError('You need to get to work on this!')

def frey_curve_model(E):
    r"""
    INPUT: An elliptic curve E defined over a field K of characteristic different from 3 or 2 having full 2-torsion.
    OUTPUT: Integeral elemetns A, B of the base field of E such that E is isomorphic to the Frey curve y^2 = x(x + A)(x - B), together with an isomorphism.
    """
    raise NotIMplementedError('Need to work on this! Figure out how to output rational maps!')

def ec_isomorphisms(E,F,JustOne = False):
    r"""
    Returns on or all isomoprhisms between two elliptic curves. Assumes that characteristic of the ground field is at least 2 or 3, but work for funciton fields, unlike the current implementation in Sage.
    """
    #from ell_generic import is_EllipticCurve
    #if not is_EllipticCurve(E) or not is_EllipticCurve(F):
        #raise ValueError, "I was promised elliptic curves. Not cool man!"
    K = E.base_ring()
    jE = E.j_invariant()
    if jE != F.j_invariant():
        if JustOne: return None
        return []
    
    from sage.rings.all import PolynomialRing
    x = PolynomialRing(QQ, 'x').gen()
    
    a1E, a2E, a3E, a4E, a6E = E.ainvs()
    a1F, a2F, a3F, a4F, a6F = F.ainvs()
    
    char = K.characteristic()
    
    if char == 2 or char == 3:
        raise NotImplementedError("Too Lazy to Worry About Positive Characteristic Just Yet!")
    c4E,c6E = E.c_invariants()
    c4F,c6F = F.c_invariants()
    
    if jE == 0:
        raise NotImplementedError("I have not worked this out for non-ridgid curves... yet!")
    elif jE == 1728:
        raise NotImplementedError("I have not worked this out for non-ridgid curves... yet!")
    else:
        m,ratio = 2,(c6E*c4F)/(c6F*c4E)
    fac = factor(ratio)
    fac_list = list(fac)
    fac_unit = QQ(fac.unit())
    #print fac_unit
    #print parent(fac_unit)
    unit_list = (x^m - QQ(fac_unit)).roots(multiplicities=False)
    sq_bool = True
    fac_rt = K(1)
    for P in fac_list:
        if not sq_bool:
            break
        if P[1] % m != 0:
            sq_bool = False
        else:
            fac_rt = fac_rt*P[0]^(ZZ(P[1]/m))
    ans = []
    if sq_bool:
        ulist = [fac_rt*u for u in unit_list]
    else:
        return None
    for u in ulist:
        s = (a1F*u - a1E)/2
        r = (a2F*u^2 + a1E*s + s^2 - a2E)/3
        t = (a3F*u^3 - a1E*r - a3E)/2
        if JustOne: return (u,r,s,t)
        ans.append((u,r,s,t))
    if JustOne: return None
    ans.sort()
    return ans


        
