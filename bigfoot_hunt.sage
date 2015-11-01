from sage.all import *

var('t1,t0')

class EllipticSurface:
    """
        A fibration of elliptic curves over P^1.
    """
    def __init__(self, ais):
        self.a1 = ais[0]
        self.a2 = ais[1]
        self.a3 = ais[2]
        self.a4 = ais[3]
        self.a6 = ais[4]

    def minimal_model(self)
        """
            Returns a model for self as an elliptic curve over \Q(P^1). In this case, each Ai has degree at most chi * i, where chi is the arithmet genus of the surface.
        """
        raise NotImplementedError("Need to adapt Tate's Algorithm")
        self.minimal_module = [A1,A2,A3,A4,A6]
        return self.minimal_model

    def arithmetic_genus(self):
        #raise NotImplementedError("Need to take the maximum")
        weights = [1,2,3,4,6]
        return max([self.minimal_model()[k].degree()/weights[k] for k in [0,1,2,3,4]])

    def TrivialLattice(self):
        raise NotImplementedError()

A_z2z8 = (2*t1*t0)^4
C_z2z8 = (t0^2 - t1^2)^4
E_z2z8 = EllipticSurface([0,A_z2z8 + A_z2z8,0,A_z2z8*C_z2z8,0])

class PartialLSeriesTree:
    """
        primes_and_lists -- A dictionary whose keys are a finite set of primes and whose values are a lists of admissible local factors for LSeries.
        
        admissible_root_numbers -- A list of admissible signs for the funcitonal equation. Set to [-1, 1] by default.

    """
    def __init__(self, primes_and_lists, admissible_root_numbers = [-1,1]):
        self.primes = primes_and_lists.keys()
        self.smoothness = max(self.primes)
        root_number_set = Set(admissible_root_numbers)
        if root_number_set == Set([1]):
            self.root_number = 1
        elif root_number_set == Set([-1]):
            self.root_number = -1
        elif root_number_set != Set([-1,1]):
            raise ValueError("admmisble_root_numbers must be [-1],[1], or [-1,1]")

    def breadth_first_itterator(self):
        raise NotImplementedError("Learn about Tree Searches")

class PartialLSeries:

    def __init__(self, dict, degree = 2):
        """
            dict -- A dictionary of the form {(p1,Ap1), ... , (pn,Apn)}
            degree -- The degree of the L-function. For now, we are sticking
                      with this, since we're recording Ap's not polynomials.
        """
        if not degree == 2:
            raise NotImplementedError("Only working with elliptic curves right now")
        self.primes = dict.keys()
        self.smoothness = max(self.primes)
        self.missing_primes = [p for p in prime_range(self.smoothness) if not p in self.primes]
        if self.missing_primes != []:
            raise NotImplementedError("Decide how to grade uninsulated data")

    def well_done(self):
        """
            Returns a slight modifictation of mean of the "well done" data defined by Mazur and Stein associated to the local data of a partial L-function for an elliptic curve. As X gets larger, this is conjectured to converge to the analytic rank of E.
        """
        X = self.smoothness
        return sum([RDF(-dict[p]*log(p)/(p*log(X))) for p in prime_range(X+1)])

    def medium_rare(self):
        raise NotImplementedError()

    def raw_bias(self):
        raise NotImplementedError()



def get_specialiazation(mbar,nbar,N):
    #Be very careful about when n = m = p, this could loop.
    """
    """
    raise NotImplementedError(" Need to prove a theorem about split vs non-split for this case before removing this error.")
    # I think Edray proved the necessary theorem back in 2008, since in this case p divides n*m.
    if gcd(m,n) == 1:
        return (m,n)
    else:
        m1 = m; n1 = n
        while (m1/n1 + 1)^2 > 2:
            m1 += p
            if gcd(m1,n1) == 1:
                return (m1,n1)
        m1 = m
        n1 += p

def trace_table(p):
    """
        Input: A prime p > 2
        
        Output: A table pairs (m mod p, n mod p) for which the curve
        E_{m,n} : Y^2 = X*(X + (2*m*n)^4)*(X+ (n^2 - m^2)^4)
        has good reduction at p together with the value of A_p.
    """
    aplist = [-1,1]
    aplist += trace_range(p,16)
    aplist.sort()
    apdict = {}
    for Ap in aplist:
        apdict[Ap] = []
    for t in range(-1,p):
        if t == -1:
            m0 = 1; n0 = 0
        else:
            m0 = t; n0 = 1
        #Fundamental domain for the action of
        #Aut(Z/2 x Z/8) is 0 < m1 < (sqrt{2} - 1) n1
        m1 = m0 + p^2
        n1 = n0 + p^3
        while (m1/n1 + 1)^2 > 2:
            n1 += p
        A = (n1^2 - m1^2)^4
        C = (2*n1*m1)^4
        E = EllipticCurve([0,A+C,0,A*C,0])
        Ap = E.ap(p)
        apdict[Ap].append((m0,n0))
    return apdict

def trace_range(p,M):
    """
        INPUT: a prime p and an intger M.
        OUTPUT: The possible values of A_p for an elliptic curve over Q with good reduction at p such that the image of E(Q)_tors mod p has order divisible by M.
    """
    # A_p = p + 1 - b*M where b >= 0
    # Start big
    # Keep going down until A_p < -2*sqrt(p)
    upper = floor(2*sqrt(p))
    starter = min(upper, p+1 - M)
    lower = -upper
    dial_down = (starter - (p + 1)) % M
    Ap = starter - dial_down
    admissible_Aps = []
    while Ap >= lower:
        admissible_Aps.append(Ap)
        Ap -= M
    return admissible_Aps
