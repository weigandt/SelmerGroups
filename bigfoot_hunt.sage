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

    def arithmetic_genus(self):
        raise NotImplementedError("Need to implement minimal model first!")

    def minimal_model(self)
        raise NotImplementedError("Need to adapt Tate's Algorithm")

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




class PartialLSeries:

    def __init__(self, primes_and_Aps):
        self.primes = primes_and_Aps.keys()
        self.smoothness = max(self.primes)

    def.insulator(self)
        return prod([p for p in prime_range(self.smoothness) if not p in self.primes])


def get_coprime(m,n,p):
    #This is a subruntine to trace_table which gets a coprime pair (m2,n2) of lowest (dictionary) height.
    raise NotImplementedError("This will loop forever when n = m = p. Need to prove a theorem about split vs non-split for this case before removing this error.")
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
    #Input: A prime p >= 11
    #Output: A table pairs (m mod p, n mod p) for which the curve
    # E_{m,n} : Y^2 = X*(X + (2*m*n)^4)*(X+ (n^2 - m^2)^4)
    # has good reduction at p together with the value of A_p.
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
