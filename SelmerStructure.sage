def fixit(f):
    R = parent(f)
    f = f.monic()
    d = f.degree()
    b = lcm([a.denominator() for a in f.coeffs()])
    f = R((b^d*f(x)).subs(x = x/b))
    return f

class Places(object):
    def __init__(self, ground_field):
        pass

class GaloisModule(object):
    def __init__(self, defining_data, ground_field):
        self.defining_data = defining_data
        self.ground_field = ground_field

    def __repn__(self):
        pass

    def CartierDual(self):
        pass

    def filtration(self, level):
        pass

    def conductor(self):
        pass

    def localization(self):
        pass

    def etale_envolope(self):
        """Takes a Galois Module M and returns a triple (L,n,w) where L is an etale algebra and w is an embedding of M into mu_n(L).
        """
        pass

class LocalGaloisModule(GaloisModule):
    def __init__(self, defining_data, ground_field, place):
        pass



class FiniteFlatGroupScheme(GaloisModule):
    pass

#This class may already be defined by Robert Miller in sage, at least in an ad hoc way.
class EtaleAlgebra(object):
    def __init__(self, defining_polynomial):
        if not defining_polynomial.is_separable():
            raise RuntimeError("Defining polynomial of an ETALE algebra must be SEPARABLE!")
        else:
            self.defining_polynomial = defining_polynomial

class GaloisCohomology(object):
    def __init__(self, M):
        self.coeff_module = M
        self.ground_field = M.ground_field

class SelmerStructure(object):
    '''Takes a Galois module M and a dictionary whose keys are places v of the ground field of M and whose values are submodules of H^1(K_v, M)
    '''
    def __init__(self, M, local_conditions):
        self.module = M
        self.local_conditions = local_conditions

    def SelmerModule(self):
        pass

    def ShafarevichTateModule(self):
        pass n n

class ClassicalSelmerStructure(SelmerStructure):
    pass

class SquareClass(object):
    def __init__(self, d):
        self.representative = d

    def ThreeFold

class QuarticPolynomial(object):
    def __init__(self, a0,a1,a2,a3,a4):
        self.a0 = a0
        self.a1 = a1
        salf.a2 = a2
        self.a3 = a3
        self.a4 = a4

    def S4_extension(self):
        R.<x> = QQ[]
        f = self.a0*x^4 + self.a1*x^3 + self.a2*x^2 + self.a3*x + self.a4
        f = fixit(f)
        F.<e> = QQ.extension(f)
        self.S4_extension = F
        self.an_S4_primitive_element = e
        return S4_extension(F)

class ThreefoldPoint(object)
    def __init__(self, E_star, P):
        self.E_star = E_star
        self.point = P
        self.E = P.curve()



