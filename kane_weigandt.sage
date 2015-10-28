r"""
################################################################################
#
#           Copyright (C) 2013 Daniel Kane, James Weigandt
#   Distributed under the terms of the GNU General Public License
#
#       This code is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#       General Pulic License for more details.
#
#   The full text of the GPL is avialable at:   http://www.gnu.org/licenses/
#
################################################################################
Version 0.1

This file contains algorithms to compute the the dimension of the 2-Selmer
groups of elliptic curves over $\QQ$ with trivial mod 2 representation. The
algorithms are based on:

Sir Peter Swinnerton-Dyer, The effect of twisting on the $2$-Selmer group,
Mathematical Proceedings of the Cambridge Philosphical Society, Volume 145,
Issue 03, November 2008, pp513-526.

and were implemented by the second author after they were patiently explained by
the first author. Any errors and ineffciencies are the complete responsiblity of
the second author.

Contact: weigandt37 AT gmail.com to report bugs, ineffciencies, or bad coding.

After you attach this file into a working session of sage, the function
AbecedarianTwoSelmerRank(A,C) will quickly return the dimension of
the 2-Selmer group S^(2)(E/Q) of the elliptic curve.

                    E(A,C) : y^2 = x(x + A)(x + C)

If you have already factored A*C*(C-A) you can feed the funciton a third
parameter S which is a list of integers consisting of -1 followed by those
distinct primes dividing A*C*(C-A).

EXAMPLE:

sage: %attach kane_weigandt.sage # Do this once while in a directory containing kane_weigandt.sage
sage: AbecedarianTwoSelmerRank(1,4)
2

TODOs:

- Automatically check (and correct) a passed set of primes.
- Improve Documentation.
- Obtain generators of S^2(E/Q) inside Q(S,2) x Q(S,2) x Q(S,2).
- Implement the Cassels pairing building off of conic functionality in sage.


"""
from sage.all import *

#Local Hilbert Codes for H^1(Q,(Z/2Z))
def local_elements(n,p):
    if p == -1:
        return [ZZ(-(sign(n)-1)/2)]
    if p > 2:
        x = n.ord(p)
        n_prime = n/(p**x)
        return [ ZZ(x % 2), ZZ(-(legendre_symbol(n_prime, p)-1)/2)]
    if p == 2:
        x = n.ord(p)
        n_prime = n/(p**x)
        if n_prime % 8 == 1:
            return [ZZ(x % 2),0,0]
        if n_prime % 8 == 3:
            return [ZZ(x % 2), 1, 0]
        if n_prime % 8 == 5:
            return [ZZ(x % 2), 0, 1]
        if n_prime % 8 == 7:
            return [ZZ(x % 2), 1, 1]

#Global Hilbert Codes for (n1,n2) in H^1(Q, (Z/2Z)x(Z/2Z) ; S)
def global_elements(n1,n2,S):
    #INPUT: Two integers n1 and n2 and a list S of places [-1, 2, p1, ...]
    #OUTPUT: Concatination of outputs of local_elements
    elements=[]
    for p in S:
        elements+= local_elements(n1,p)
        elements+= local_elements(n2,p)
    return elements

#Global Hilbert Code basis for H^1(Q, (Z/2Z)x(Z/2Z) ; S)
def list_of_global_elements(S):
    #INPUT: List of p
    #OUTPUT:
    element_list=[]
    for p in S:
        element_list.append(global_elements(p,1,S))
        element_list.append(global_elements(1,p,S))
    return element_list

#Linear Algebra
def add_if_indep(L,v):
    #INPUT: A list L of F_2 vectors (as lists) and another F_2 vector (as a list)
    #OUTPUT: L if v is in the span of L, L + v otherwise
    V = matrix(GF(2),[v])[0]
    if V == 0:
        return L
    if V in matrix(GF(2),L).row_space():
        return L
    return L + [v]

#Finds least quadratic non-residue.
#TODO: Make this unnecessary.
def quad_non_res(p):
    alpha = 1
    while legendre_symbol(alpha,p) != -1:
        alpha += 1
        if alpha > p:
            raise RuntimeError("No quadratic non-residue found, this should not happen!")
    return alpha

#Local Hilbert Codes for E(Q_p)/2E(Q_p) when E : y^2 = x(x + A)(x + C)
def local_mw_basis(A,C,p):
    #INPUT: (A,C,p)
    #OUTPUT: Generators of E(Q_p)/2E(Q_p) as a d x 2d matrix over F_2
    mw_basis=[]
    mw_basis = add_if_indep(mw_basis,local_elements(A*C,p)+local_elements(A,p))
    #print mw_basis
    mw_basis = add_if_indep(mw_basis,local_elements(-A,p) + local_elements(-A*(C-A),p))
    mw_basis = add_if_indep(mw_basis,local_elements(-C,p) + local_elements(-(C-A),p))
    d = ZZ(len(local_elements(1,p)))
    #print d
    mw_dim = matrix(GF(2),mw_basis).rank()
    if p > 2:
        if mw_dim < d:
            B = C - A
            if gcd(A,C) % p == 0:
                raise RuntimeError("This should never happen, two torsion should determine local images at " + str(p))
            alpha = quad_non_res(p)
            if ZZ(A*B*C) % p != 0:
                raise RuntimeWarning("Why are you checking local images at " + str(p) + "?")
                mw_basis = add_if_indep(mw_basis, local_elements(alpha,p) + local_elements(alpha,p))
                mw_basis = add_if_indep(mw_basis, local_elements(alpha,p) + local_elements(1,p))
                mw_dim = len(mw_basis)
            if ZZ(A) % p == 0:
                mw_basis = add_if_indep(mw_basis, local_elements(alpha,p) + local_elements(alpha,p))
                mw_dim = len(mw_basis)
                if mw_dim < d:
                    if ZZ(A).ord(p) % 2 != 0:
                        raise RuntimeError("2 torsion plus (a,a,1) should have determined the local image at " + str(p) + "!")
                    if QQ(C).is_padic_square(p):
                        mw_basis = add_if_indep(mw_basis, local_elements(p,p) + local_elements(p,p))
                        mw_dim = len(mw_basis)
                        if mw_dim < d:
                            raise RuntimeError("Something has gone wrong with case 2b1 at the prime " +str(p))
                    else:
                        mw_basis = add_if_indep(mw_basis, local_elements(alpha,p) + local_elements(1,p))
                        mw_dim = len(mw_basis)
                        if mw_dim < d:
                            raise RuntimeError("Something has gone wrong with case 2b2 at the prime " + str(p))
            if ZZ(B) % p == 0:
                mw_basis = add_if_indep(mw_basis, local_elements(1,p) + local_elements(alpha,p))
                mw_dim = len(mw_basis)
                if mw_dim < d:
                    if ZZ(B).ord(p) % 2 != 0:
                        raise RuntimeError("E[2] plus (1,a,a) should have determined the local image at " +str(p) + " dividing B")
                    if QQ(-A).is_padic_square(p):
                        mw_basis = add_if_indep(mw_basis, local_elements(1,p) + local_elements(p,p))
                        mw_dim = len(mw_basis)
                        if mw_dim < d:
                            raise RuntimeError("Something has gone wrong with case 2b1 at the prime " + str(p) + " dividing B!")
                    else:
                        mw_basis = add_if_indep(mw_basis, local_elements(alpha,p) + local_elements(alpha,p))
                        mw_dim = len(mw_basis)
                        if mw_dim < d:
                            raise RuntimeError("Something has gone wrong with case 2b2 at the prime " + str(p) + " dividing B!")
            if ZZ(-C) % p == 0:
                mw_basis = add_if_indep(mw_basis, local_elements(alpha,p) + local_elements(1,p))
                mw_dim = len(mw_basis)
                if mw_dim < d:
                    if ZZ(-C).ord(p) % 2 != 0:
                        raise RuntimeError("E[2] plus (a,1,a) should have determined the local image at " + str(p) + " dividing C!")
                    if QQ(-B).is_padic_square(p):
                        mw_basis = add_if_indep(mw_basis, local_elements(p,p) + local_elements(1,p))
                        mw_dim = len(mw_basis)
                        if mw_dim < d:
                            raise RuntimeError("Something has gone wrong with the case 2b1 at the prime " + str(p) + " dividing C!")
                    else:
                        mw_basis = add_if_indep(mw_basis, local_elements(1,p) + local_elements(alpha,p))
                        mw_dim = len(mw_basis)
                        if mw_dim < d:
                            raise RuntimeError("Something has gone wrong with the case 2b2 at the prime " + str(p) + " dividing C!")
    N=1
    while N*(N+A)*(N+C)*(N-A)*(N-C)==0:
        N=N+1/4
    while mw_dim < d:
        v1 = local_elements(ZZ(4*N),p)
        v2 = local_elements(ZZ(4*(N+A)),p)
        v3 = local_elements(ZZ(4*(N+C)),p)
        if matrix(GF(2),v1) + matrix(GF(2),v2) + matrix(GF(2),v3) == 0:
            mw_basis = add_if_indep(mw_basis,v1 + v2)
            mw_dim = len(mw_basis)
        N = N + 1/4
        while N*(N+A)*(N+C)*(N-A)*(N-C)==0:
            N=N+1/4
    #print N
    return mw_basis

#Creates the big matrix whose kernel is isomophic to Sel(2)(E/\QQ)
def get_underlying_matrix(A,C,S=None):
    #INPUT: A,C integers defining our Frey Curve y^2 = x(x+A)(x+C), 
    #        S = [-1, 2, p1, p2, ...] list of "primes"
    #OUTPUT: A Big F_2 matrix for selmer group computation
    if S is None:
        S = [-1] + ZZ(A*C*(A-C)).support()
    M = list_of_global_elements(S)
    big_dim = len(M[0])
    leading_zeros=[]
    for p in S:
        L=local_mw_basis(A,C,p)
        for v_small in L:
            v_big = leading_zeros + v_small
            while len(v_big) < big_dim:
                v_big.append(0)
            M.append(v_big)
        #print v_big
        if p == -1:
            leading_zeros = leading_zeros + [0,0]
        if p == 2:
            leading_zeros = leading_zeros + [0,0,0,0,0,0]
        if p > 2:
            leading_zeros = leading_zeros + [0,0,0,0]
    return M

def AbecedarianTwoSelmerRank(A,C,S = None, verbose = False):
    #INPUT: A,C integers, list of places containing [-1, 2] and bad primes for the associated Frey Curve
    #OUTPUT: The 2-Selmer Rank of E : y^2 = x(x+A)(x+C).
    M=get_underlying_matrix(A,C,S)
    M=matrix(GF(2),M)
    if verbose:
        V = M.kernel()
        print "Hilbert Code and Local Image Matrix:"
        print M
        print "Kernel Representation of Selmer Group:"
        print V
    s = M.nullity()
    if verbose:
        print "Upper bound on Mordell-Weil rank is " + str(s) + "."
    return s
