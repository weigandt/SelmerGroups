from psage.ellcurve.xxx.rankbound import xxx_rankbound
import numpy as np
import numpy.linalg

bad_primes = [2,3,5,7,11,13,17,19,48463,20650099, 315574902691581877528345013999136728634663121, 376018840263193489397987439236873583997122096511452343225772113000611087671413]

def triple_xxx_rank(E,Delta,d,verbose = False, bad_primes = bad_primes):
    deltas = srange(1,Delta,d)
    x1 = np.array([1/t for t in deltas])
    y1 = np.array([xxx_rankbound(E,t,bad_primes = bad_primes) for t in deltas])
    n = len(x1)
    A = np.array([[x1[j],1] for j in range(n)])
    B = y1
    X = np.linalg.lstsq(A,B)
    if verbose:
        print X[0]
        print X[1]
    return X[0][1]

def triple_xxx_trim(Delta, d, r, input_file, progress_field, output_file):
    lines = input_file.readlines()
    curves = []
    for T in lines:
        ainvars = sage_eval(T)
        E = EllipticCurve(ainvars)
        curves.append(E)
    for E in curves:
        print("Chekcing Curve: " + "\n" + str(E.a_invariants()) + "\n" + "with D = " + str(Delta) + " and d = " + str(d))
        rg = triple_xxx_rank(E,Delta,d,verbose = True)
        if rg > r - 2:
            output_file.write(str(E.a_invariants()) + '\n')
            output_file.write(str([rg,Delta,d]) + '\n')

def xxx_trim(delta,r, input_file, progress_file, output_file):
    lines = input_file.readlines()
    curves = []
    for T in lines:
        ainvars = sage_eval(T)
        E = EllipticCurve(ainvars)
        curves.append(E)
    for E in curves:
        print("Checking Curve: " + "\n" + str(E.a_invariants()) + '\n' + "with delta = " + str(delta) +'\n')
        rb = xxx_rankbound(E,delta)
        print rb
        if rb > r - 0.001:
            output_file.write(str(E.a_invariants()) + '\n')

