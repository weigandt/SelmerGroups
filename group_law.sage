from sage.all import *

def inverse_of_point(E,P0):
    if len(P0) > 2:
        if not P0 == (0,1,0):
            raise NotImplementedError("Only works for affine points and (0,1,0)!")
        else:
            return (0,1,0)
    else:
        x0 = P0[0]
        y0 = P0[1]
        return (x0, -y0 - E.a1()*x0 - E.a3())

def add_points(E,P1,P2):
    #Computes the sum P1 + P2 on E
    if len(P1) > 2:
        if not P1 == (0,1,0):
            raise NotImplementedError("Only adds affine points and (0,1,0)!")
        else:
            if len(P2) > 2:
                if not P2 == (0,1,0):
                        raise NotImplementedError("Only adds affine points, and (0,1,0)!")
                else: 
                    return (0,1,0)
            else:
                return P2
    elif len(P2) > 2:
        if not P2 == (0,1,0):
            raise NotImplementedError("Only adds affine points and (0,1,0)!")
        else:
            return P1
    else:
        x1 = P1[0]
        y1 = P1[1]
        x2 = P2[0]
        y2 = P2[1]
        a1,a2,a3,a4,a6 = E.a_invariants()
        if (y1^2 + a1*x1*y1 + a3*y1) != (x1^3 + a2*x1^2 + a4*x1 + a6):
            raise RuntimeError("Point (" + str(x1) + ", " + str(y1) + ") is does not lie on " + str(E) + ",")
        elif (y2^2 + a1*x2*y2 + a3*y2) != (x2^3 + a2*x2^2 + a4*x2 + a6):
            raise RuntimeError("Point (" + str(x2) + ", " + str(y2) + ") is does not lie on " + str(E) + ",")
        elif x1 == x2:
            if (y1 + y2 + a1*x2 + a3) == 0:
                return (0,1,0)
            else:
                lambdan = (3*x1^2 + 2*a2*x1 + a4 - a1*y1)
                lambdad = (2*y1 + a1*x1 + a3)
                nun = (-x1^3 + a4*x1 + 2*a6 - a3*y1)
                nud = (2*y1 + a1*x1 + a3)
                lam = lambdan/lambdad
                nu = nun/nud
        else:
            lam = (y2 - y1)/(x2 - x1)
            nu = (y1*x2 - y2*x1)/(x2 - x1)
        x3 = lam^2 + a1*lam - a2 - x1 - x2
        y3 = -(lam + a1)*x3 - nu - a3
        return (x3,y3)