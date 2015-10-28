from sage.all import *

#The BinaryQuadraticForm code seems difficult to work with in sage, I will begin with my own code an integrate as necessary.

def I_J(quartic):
    a,b,c,d,e = quartic
    I = 12*a*e - 3*b*d + c^2
    J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3
    return I, J

def bqf_jac(quartic):
    I,J = I_J(quartic)
    return EllipticCurve([-27*I,-27*J])

def two_cover(quartic):
    if len(quartic) != 5:
        raise RuntimeError("Quartic " + str(quartic) + "does not represent a binary quartic form")
    I,J = I_J(quartic)
    a,b,c,d,e = quartic
    R.<X,Y,Z> = PolynomialRing(parent(a*b*c*d*e).fraction_field(), 3)
    g4 = (3*b^2 - 8*a*c)*X^4 + 4*(b*c - 6*a*d)*X^3*Y + 2*(2*c^2 - 24*a*e - 3*b*d)*X^2*Y^2 + 4*(c*d - 6*b*e)*X*Y^3 + (3*d^2 - 8*c*e)*Y^4
    g6 = (b^3 + 8*a^2*d - 4*a*b*c)*X^6 + 2*(16*a^2*e + 2*a*b*d - 4*a*c^2 + b^2*c)*X^5*Y + 5*(8*a*b*e + b^2*d - 4*a*c*d)*X^4*Y^2 + 20*(b^2*e - a*d^2)*X^3*Y^3 - 5*(8*a*d*e + b*d^2 - 4*b*c*e)*X^2*Y^4 - 2*(16*a*e^2 + 2*b*d*e - 4*c^2*e + c*d^2)*X*Y^5 - (d^3 + 8*b*e^2 - 4*c*d*e)*Y^6
    return (3*g4.subs(Y = 1)/(4*Y^2), 27*g6.subs(Y = 1)/(8*Y^3))

def three_cover(M):
    raise NotImplementedError

def four_cover(A,B):
    raise NotImplementedError

