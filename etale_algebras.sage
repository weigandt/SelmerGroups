from sage.all import *

def hyperelliptic_J2_algebra(f):
    return tuple([g[0] for g in list(f.factor())])

def etale_algebra2(E):
    #INPUT: An elliptic curve E over QQ
    #OUTPUT: A tuple of polynomials giving maps from E to H^1(QQ,E[2])
    f = E.minimal_model().short_weierstrass_model().division_polynomial(2).monic()
    return hyperelliptic_J2_algebra(f)

def etale_selmer2(f,S):
    #INPUT: An etale algebra A written as a list of irreducible polynomials, and a finite set S of primes
    #OUTPUT: A basis for the (2,S)-Selmer group of this etale algebra.
    if parent(f) = sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field_with_category:
        f = f.minimal_model().short_weierstrass_model().division_polynomial(2).monic()
    A = hyperelliptic_J2_algebra(f)
    bases = []
    for g in A:
        if g.degree() == 1:
            gbasis = [QQ(-1)]
            for p in S:
                gbasis.append(QQ(p))
            bases.append(gbasis)
        else:
            L.<e> = g.root_field()
            SL = L.primes_above(prod(S))
            bases.append(L.selmer_group(SL,2))
    big_basis = []
    for it in range(len(bases)):
        for B in bases[it]:
            element = [QQ(1) for g in A]
            element[it] = B
            element = tuple(element)
            big_basis.append(element)
    return big_basis

#Write and Etale Algebra Class and overload * for this.
def etale_product(dlist,nofields = 1):
    if len(dlist) == 0:
        return tuple([QQ(1) for i in range(nofields)])
    return tuple([prod([d[i] for d in dlist]) for i in range(len(dlist[0]))])

#Find a way to use linear algebra for this. The norm should be seen as a linear map.
def norm_kernel(basis, nofields = 1):
    ker = []
    for dlist in subsets(basis):
        D = etale_product(dlist, nofields)
        if QQ(prod([QQ(d.norm()) for d in D])).is_square():
            ker.append(D)
    return ker

def ker_of_norm(E):
    m = 2
    S_ground = ZZ(m*E.conductor()).support()
    L, e, f = two_div_field(E)
    S = L.primes_above(S_ground)
    unr_gens = L.selmer_group(S,m)
    ker = []
    for d in subsets(unr_gens): # Specific to m = 2
        D = L(prod(d))
        if QQ(D.norm()).is_square():
            ker.append(D)
    return ker, L, e
