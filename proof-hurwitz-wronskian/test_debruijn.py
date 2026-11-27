import numpy as np
import sympy as sp

def wronskian_roots(rootsP, rootsQ):
    t=sp.symbols('t')
    P=1
    for r in rootsP: P*=(t-r)
    Q=1
    for r in rootsQ: Q*=(t-r)
    W=sp.expand(P*sp.diff(Q,t)-sp.diff(P,t)*Q)
    poly=sp.Poly(W,t)
    deg=poly.degree()
    if deg<0: return []
    cs=[int(poly.nth(i)) for i in range(deg+1)]
    desc=cs[::-1]
    return np.roots(desc)

# interlacing roots: P roots -1, -3? Let's test random
import random
random.seed(0)
for trial in range(5):
    # generate sorted negative roots for P degree 3 and Q degree 4 interlacing
    pts=sorted(-random.random()*10 for _ in range(7))
    # P roots: pts[1], pts[3], pts[5]; Q roots: pts[0],pts[2],pts[4],pts[6]
    P=[pts[1],pts[3],pts[5]]
    Q=[pts[0],pts[2],pts[4],pts[6]]
    rts=wronskian_roots(P,Q)
    maxre=max(r.real for r in rts)
    print("trial",trial,"max real part", maxre, "all neg?", all(r.real<1e-9 for r in rts))
