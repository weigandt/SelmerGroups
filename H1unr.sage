from sage.all import *

##########################################################################
# Jamie: This is from my 2-descent code. It's not as fast as it could be
#         but it should get the job done for what you need it for.
#         It also works regardless of the mod 2 rep'n
##########################################################################


#################################################
#
#           Basic 2-division field stuff
#
#################################################

def two_div_field(E):
    a4 = E.short_weierstrass_model().a4()
    a6 = E.short_weierstrass_model().a6()
    f = x^3 + a4*x + a6
    L.<alpha> = QQ.extension(f)
    return L, alpha, f

def fixit(f):
    R = parent(f)
    f = f.monic()
    d = f.degree()
    b = lcm([a.denominator() for a in f.coeffs()])
    f = R((b^d*f(x)).subs(x = x/b))
    return f

######################################################################
#
# Computes Unramified Cohomology Classes
#
######################################################################

# In general, if f is any separable polynomial of odd degree, then the hyperelliptic curve C : y^2 = f(x) has
#  H^1(Q, Jac(C)[2]) = ker [ A^\times / (A^\times)^2 ---> Q^\times/(Q^\times)^2 ]
#    where A = Q[x]/(f(x)) is the so-called etale-algebra associated to Jac(C)[2].
#  We also have the "Selmer Group", A(S,2) which sits inside the exact sequence
#  0 --> O_{A,S}^\times/ (O_{A,S}^\times)^2 --> A(S,2) --> Cl(O_{A,S})[2] --> 0
# Finally, we compute H^1(Q,Jac(C)[2];S) = H^1(Q,Jac(C)[2]) \cap A(S,2)
#NOTE: There should be a faster way to do these calculations by using standard linear algebra algorithms for
#      computing the kernel of a linear map, and the intersection of two subsapces found in Cohen's Book.

def hyperelliptic_J2_algebra(f):
    #TODO: Include some verification that f is separable, and possibly of odd degree.
    return tuple([g[0] for g in list(f.factor())])

def etalg(E):
    #INPUT: An elliptic curve E over QQ
    E = E.minimal_model()
    var('x')
    R = QQ[x]
    f = R(x^3 -27*E.c4()*x - 54*E.c6())
    fac = list(f.factor())
    etalg = [g[0] for g in fac]
    return tuple(etalg)

def etsel(A,S):
    #INPUT: An etale algebra A written as a list of irreducible polynomials, and a finite set S of primes
    #OUTPUT: A basis for the (2,S)-Selmer group of this etale algebra.
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

#Write a square class class an overload multiplication with this.
def etprod(dlist,nofields = 1):

    if len(dlist) == 0:
        return tuple([QQ(1) for i in range(nofields)])
    return tuple([prod([d[i] for d in dlist]) for i in range(len(dlist[0]))])

#This is slow. It calls product WAY too many times. There should be a more direct way to compute things here.
#Consider using the Snake Lemma with norm maps.
def norm_kernel(basis, nofields = 1):
    ker = []
    for dlist in subsets(basis):
        D = etprod(dlist, nofields)
        if QQ(prod([QQ(d.norm()) for d in D])).is_square():
            ker.append(D)
    return ker

#### Returns list of elements of H^1(Q,E[2];S) as tuples of number field elements.

def unramified_square_classes(E): # I think, Mazur-Rubin call this "unrestricted at S Selmer Module"
    A = etalg(E)
    S = ZZ(2*E.conductor()).support()
    Sel2A = etsel(A,S)
    t = E.two_torsion_rank()
    return norm_kernel(Sel2A, nofields = t + 1)

def UnrestrictedSelmerGroup(E):
    """
    INPUT:
    OUTPUT:
    """
    return unramified_square_classes(E)

def associated_quartic(d):
    """
    INPUT: d = (d_i)_{i = 1}^{m} where d_i is in some number field k_i.
    OUTPUT: Quartic polynomials associated to d via different interprestations of Selmer groups
    """
    if len(d) != 1:
        raise NotImplementedError("Stop being a coward, Jamie.")
    d0 = d[0]
    k = parent(d0)
    K.<eK> = k.galois_closure()
    em = k.embeddings(K)
    d1 = em[0](d0)
    d2 = em[1](d0)
    RK.<x> = K[]
    F.<eF> = K.extension((x^2 + RK(d1) - RK(d2))^2 - 4*x^2*RK(d1))
    Fabs.<eFabs> = F.absolute_field()
    Ls = [data[0] for data in Fabs.subfields(degree = 4)]
    return Ls[0]

def EdraysSection(d):
    if len(d) != 1:
            raise NotImplementedError("Think about this section for non-surjective mod 2's")
    d0 = d[0]
    tr = d0.trace()
    n = d0.norm()
    a0 = tr^2 -4*n*(1/d0).trace()
    a1 = -8*sqrt(n)
    a2 = -2*tr
    a3 = 0
    a4 = 1
    I = a2^2 - 3*a1*a3 + 12*a0*a4
    J = -2*a2^3 + 9*a1*a2*a3 - 27*a0*a3^2 - 27*a1^2*a4 + 72*a0*a2*a4
    E = EllipticCurve([-I/3,-J/27])
    P = E(-2*a2/3, 3*a1/3)
    F.<eF> = NumberField(x^4 + a3*x^3 + a2*x^2 + a1*x + a0)
    return P, (a0,a1,a2,a3,a4), F, E, I, J 

def classgroup_rankbound(E):
    A = etalg(E)
    S = ZZ(2*E.conductor()).support()
    t = E.two_torsion_rank()
    Sel2A = etsel(A,S)
    H1unr = norm_kernel(Sel2A, nofields = t + 1)
    num = len(H1unr)
    r0 = ZZ(num).ord(2)
    if num != 2^r0:
        raise RuntimeError("Number of cohomology classes should be a power of 2!")
    return r0 - t

def test_hypothesisA(E):
    if E.two_torsion_rank() != 0:
        return True
    L,e,f = two_div_field(E)
    if L.galois_closure('eLgc').degree() != 6:
        return True
        #Really should be not implemented error, but the hypothesis isn't implemented either.
    ds = [d for d in UnrestrictedSelmerGroup(E) if d[0] != 1]
    Sold = ZZ(2*E.conductor()).support()
    for d in ds:
        Snew = ZZ(EdraysSection(d).curve().conductor()).support()
        if union(Snew,Sold) != Sold:
            print E.cremona_label(), d, Sold, Snew
            return False
    return True

def test_hypothesisB(E):
    if E.two_torsion_rank() != 0:
        return True
    L,e,f = two_div_field(E)
    if L.galois_closure('eLgc').degree() != 6:
        return True
        #Really should be not implemented error, but the hypothesis isn't implemented either.
    ds = [d for d in UnrestrictedSelmerGroup(E) if d[0] != 1]
    Sold = Set(ZZ(2*E.conductor()).support())
    for d in ds:
        a0,a1,a2,a3,a4 = EdraysSection(d)[1]
        g = a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
        F.<eF> = NumberField(g)
        Snew = Set(F.discriminant().support())
        if not Snew.issubset(Sold):
            print E.cremona_label(), d, Sold, Snew
            return False
    return True


def test_hypothesisC(E):
    if E.two_torsion_rank()!=0:
        return True
    L,e,f=two_div_field(E)
    if L.galois_closure('eLgc').degree() != 6:
        return True
    E = E.short_weierstrass_model()
    A = E.a4()
    B = E.a6()
    Delta = -16*(4*A^3 + 27*B^2)
    x = PolynomialRing(RationalField(), 'x').gen()
    f = x^4 - 4*Delta*x - 12*A*Delta
    if f.is_irreducible():
        return True
    dss=[d for d in UnrestrictedSelmerGroup(E) if d[0] != 1]
    for d in dss:
        F1 = associated_quartic(d)
        F2 = BFextn(EdraysSection(d)[0])
        F4 = EdraysSection(d)[2]
        if not F1.is_isomorphic(F4) or not F1.is_isomorphic(F2):
            print E.cremona_label(),d
            return False
        return True
    
def test_hypothesisD(E):
    if E.two_torsion_rank()!=0:
        return False
    L,e,f=two_div_field(E)
    if L.galois_closure('eLgc').degree() != 6:
        return False
    E = E.short_weierstrass_model()
    A = E.a4()
    B = E.a6()
    Delta = -16*(4*A^3 + 27*B^2)
    x = PolynomialRing(RationalField(), 'x').gen()
    f = x^4 - 4*Delta*x - 12*A*Delta
    if f.is_irreducible():
        return True
    else:
        return False

def show_isom(E):
    dss=[d for d in UnrestrictedSelmerGroup(E) if d[0] != 1]
    for d in dss:
        d = d[0]
        t = d.trace()
        n = d.norm()
        g = (1/d).trace()
        Q = x^4 - 32*g*n*x^2 - 512*n*x + 256*g^2*n^2 - 1024*n*t
        q = x^4-2*t*x^2-ZZ(8*sqrt(n))*x+(t^2-4*n*g)
        Q = QQ[x](Q)
        q = QQ[x](q)
        K.<A> = NumberField(Q)
        k.<a> = NumberField(q)
        fK = K.pari_polynomial()
        fk = k.pari_polynomial()
        print t,n,g, fK.nfisisom(fk), fk.nfisisom(fK)
        
        
def S4_quartic_field(d):
    a0,a1,a2,a3,a4 = EdraysSection(d)[1]
    g = a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
    F.<eF> = NumberField(g)
    return F

def TateModuleS4(E):
    E = E.short_weierstrass_model()
    A = E.a4()
    B = E.a6()
    Delta = -16*(4*A^3 + 27*B^2)
    f = x^4 - 4*Delta*x - 12*A*Delta
    F.<eF> = NumberField(f)
    return F

def BFextn(P):
    xP = P[0]/P[2]
    a1,a2,a3,a4,a6 = P.curve().a_invariants()
    F.<eF> = NumberField(x^4 - 4*xP*x^3 -(2*a4+a1*a3+a1^2*xP + 4*a2*xP)*x^2 - 2*(a3^2 + 4*a6 + 2*a4*xP + a1*a3*xP)*x - (a1^2*a6 +4*a2*a6 - a1*a3*a4 + a2*a3^2 -a4^2 + a3^2*xP + 4*a6*xP))
    return F


                         