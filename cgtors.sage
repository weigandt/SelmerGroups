from sage.all import *

def store_factor_base(T,OK):
    L = OK.number_field()
    B = {}
    prime_list = []
    for p in prime_range(T):
        B[p] = []
        for P in L.primes_above(p):
            if P.norm() < T:
                B[p].append(P)
                prime_list.append(P)
    return (OK,T,B,prime_list)

#Returning a prime_list is important because it will allow us to order matrix entries by residue characteristic. 
##Used to extend a factor base, back when the 4th peice of data was a list of
## primes with small characteristic and large norm. Never seemed useful, but
## might be in a faster calculation of Belabas' bound.
#def extend_factor_base(B,T):
#    OK = B[0]
#    L = OK.number_field()
#    lilT = B[1]
#    lilB = B[2]
#    lilBextra = B[3]
#    bigB = lilB
#    bigBextra = {}
#    for p in prime_range(lilT):
#        bigBextra[p] = []
#        for P in lilBextra[p]:
#            if P.norm() < T:
#                bigB[p].append(P)
#            else:
#                bigBextra.append(P)
#    for p in prime_range(lilT,T):
#        bigB[p] = []
#        bigBextra[p] = []
#        for P in L.primes_above(p):
#            if P.norm() < T:
#                bigB[p].append(P)
#            else:
#                bigBextra.append(P)
#    return (OK, T, bigB, bigBextra)

def BDF_sum(T,B, Bound):
    L = B[0].number_field()
    c1 = pi^2/2
    c2 = 4*catalan
    n = L.degree()
    r1 = len(L.embeddings(RR))
    S0 = -(n*c1 + r1*c2)/log(T)
    S = float(S0)
    for p in prime_range(T):
        for P in B[2][p]:
            Np = P.norm()
            m = 1
            while Np^m < T:
                S += float(2*log(Np)/Np^(m/2) * ( 1- log(Np^m)/log(T)))
                m += 1
        if S > Bound:
            return True
        if floor(Li(p)) % 20 == 0:
            print S
    return false

def Bach_bound(OK):
    L = OK.number_field()
    Del = L.disc(OK.basis())
    Tbach = 12*log(Del.abs())^2
    return Tbach

def BDF_bound(OK):
    #INPUT:     The maximal order OK of a number field L
    #OUTPUT:    A bound T such that the ideal class group of OK is generated by prime
    #           ideals of norm at most T assuming the GRH for the Dedekind zeta
    #           function of L.
    L = OK.number_field()
    r1 = len(L.embeddings(RR))
    Del = L.disc(OK.basis())
    n = L.degree()
    D = log(Del.abs()) - n*(euler_gamma + log(8*pi)) - r1*pi/2
    print D
    print RR(D)
    Tbach = Bach_bound(OK)
    B = store_factor_base(Tbach,OK) #May take a while, 100 sec for L28
    Tgood = ceil(Tbach)
    Tbad = 0
    while Tgood - Tbad > 2:
        Ttry = floor(Tgood + Tbad)/2
        print "Trying " + str(Ttry)
        check = BDF_sum(Ttry, B, D)
        if check:
            Tgood = ceil(Ttry)
        else:
            Tbad = floor(Ttry)
    return Tgood

def prime_list(bigB):
    #INPUT: A factor base (OK,T,B)
    #OUTPUT: A list of all the prime ideals of OK in B
    return [P for P in bigB[2].itervalues()]

#Slow compared to factoring the ideal, but works well if factorization is slow, returns aux. ideal.
def record_valuations_mod_ell(alpha,B,ell = 2):
    #INPUT: an element alpha of a number field, a factor base (OK,T,B), and a prime ell
    #OUPUT: a list of valuations of alpha at primes P in Plist mod ell
    #       and the norm of the ideal left over
    primes_to_hash = []
    alphavec = []
    I = alpha*B[0]
    NI = I.norm()
    for p in B[2].iterkeys():
        if NI.ord(p) != 0:
            for P in B[2][p]:
                ordP = I.valuation(P)
                alphavec.append(GF(ell)(ordP))
                I = I*P^(-ordP)
        else:
            for P in B[2][p]:
                alphavec.append(GF(ell)(0))
    return alphavec, I

#Fast if the ideal can be factored quickly. returnse set of primes
def Bproj_mod_ell(alpha,B,ell=2,rnum = 0):
    #INPUT: - an element alpha to consider the princpal ideal alpha*OK
    #       - A factor base (O,T,dict)
    #       - An rational prime ell
    #OUTPUT: - The dictionary describing the image of (alpha) in Fp^B
    #           for initialization of a sparse matrix
    #        - A list of primes P outside for which ord_P(alpha) is odd.
    I = alpha*B[0]
    Plist = B[3]
    NI = I.norm()
    alphaitems = []
    primes_to_hash = []
    for p in NI.support():
        if p < B[1]:
            for P in B[2][p]:
                ordP = I.valuation(P)
                alphaitems.append(((rnum,Plist.index(P)),GF(ell)(ordP)))
        if I.norm() % p == 0:
            for P in B[0].number_field().primes_above(p):
                if P.norm() > B[1]:
                    ordP = I.valuation(P)
                    if GF(ell)(ordP) != 0:
                        primes_to_hash.append(P)
    return alphaitems, primes_to_hash

def trivial_relations(B,ell=2,verbose = False):
    #INPUT: a factor base (OK,T, dict) of primes of OK up to norm T, a prime ell
    #OUTPUT: a list of trivial relations mod ell form considering p*OK
    ncols = sum([len(B[2][p]) for p in B[2].iterkeys()])
    V = span(Matrix(GF(2), {},nrows = 1, ncols = ncols))
    for p in B[2].iterkeys():
        Bproj_data = Bproj_mod_ell(p,B,ell)
        if Bproj_data[1] == []:
            Vp = span(Matrix(GF(2), Bproj_data[0], nrows = 1, ncols = ncols))
            V = V + Vp
            if verbose:
                print p,V
    return V

def trivial_dictionary(B,ell = 2):
    #INPUT: - factor base B = (OK,T,dict, list)
    #       - a rational prime ell
    #OUTPUT: -a dictionary describing the trivial realtions mod ell coming from the the ideals p*O for p < T
    Tdict = {}
    rownum = 0
    for p in prime_range(B[1]):
        pitems, bad = Bproj_mod_ell(p,B,ell,rnum = rownum)
        if bad == []:
            for it in pitems:
                Tdict[it[0]] = it[1]
            print p, rownum
            rownum += 1
    print("Trivial upper bound on GF(" + str(ell) + ") rank is " + str(len(B[3]) - rownum) + ".")
    return Tdict
