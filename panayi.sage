from sage.all import *

def unif(pnf):
    return pnf.gens()[-1]

def panayi_sharp(pnf,f):
    pi = unif(pnf)
    vpf = min([pnf.number_field()(c).ord(p) for c in f.coeffs()])
    #print vpf,f, f.coeffs()
    return f/pi^vpf

def residual_roots(pnf,c):
    K = pnf.number_field()
    x = K.polynomial_ring().gens()[0]
    roots = []
    for beta in pnf.residues():
        if c.subs(x = beta) in pnf:
            roots.append(K(beta))
    return roots
    
def reduction_degree(pnf,h):
    h = panayi_sharp(pnf,h)
    cfs = h.coeffs()
    red_data = [(cfs[i].ord(pnf),i) for i in range(h.degree() + 1)]
    d = 0
    for D in red_data:
        if D[0] < 0:
            raise RuntimeError("Sharp function is off! h-sharp should not be " + str(h))
        elif D[0] == 0:
            if D[1] > d:
                d = D[1]
    return d

def replacements(pnf,c):
    K = pnf.number_field()
    x = K.polynomial_ring().gens()[0]
    R = residual_roots(pnf,c)
    pi = unif(pnf)
    replist = []
    m_incr = 0
    for beta in R:
        h = c.subs(x = pi*x + beta)
        h = panayi_sharp(pnf,h)
        d = reduction_degree(pnf,h)
        if d == 1:
            m_incr += 1
            print ("Root found: " + str(h))
        if d > 1:
            replist.append(h)
    return replist, m_incr

def panayi_root_count(pnf,f):
    #INPUT: A prime pnf in a number field L, and a squarefree polynomial f in QQ[x].
    #OUTPUT: The number of roots of f in L_p.
    print("------------------------------")
    print("Finding roots for " + str(f) + ".")
    print("------------------------------")
    K = pnf.number_field()
    f = K.polynomial_ring()(f)
    if not f.is_squarefree():
        raise RuntimeError("Panayi Algorithm Implemented Only for Squarefree f")
    pi = unif(pnf)
    F = pnf.residue_field()
    Rp = F.polynomial_ring()    
    C = [panayi_sharp(pnf,f)]
    print("Initial Factors: " + str(C))
    m = 0
    print("m = 0")
    while len(C) > 0:
        Ctemp = []
        for c in C:
            rep, minc = replacements(pnf,c)
            Ctemp = Ctemp + rep
            m = m + minc
        C = Ctemp
        print("New Factors: " + str(C))
        print("m = " +str(m))
    return m

def check_2adic(L,prec,gens):
    p = L.primes_above([2])[0]
    I = p^prec
    output = []
    for B in subsets(gens):
        d = prod(B)
        output.append(sum([(T^2 - d in I) for T in I.residues()]))
    return output


