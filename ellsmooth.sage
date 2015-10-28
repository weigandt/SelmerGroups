from sage.all import *

def is_primitive(E):
    F = E.minimal_quadratic_twist()[0]
    return abs(E.minimal_model().discriminant()) == abs(F.minimal_model().discriminant())

def kap(E):
    S = max(E.conductor().support())
    Del = abs(E.minimal_model().discriminant())
    kp = n(log(S)/log(log(Del)),12)
    return kp

def smoplt(X):
    data = []
    for E in cremona_curves(1..X):
        if is_primitive(E):
            N = E.conductor()
            kp = kap(E)
            if kp < 1/2:
                print N, kp
            data.append([N,kp])
    return list_plot(data, size = 1, dpi = 100)

