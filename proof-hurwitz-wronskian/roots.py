import sys, math
sys.path.insert(0, 'C:/Users/Abdelli/Documents/ordered-big-blocks')
from scan import ordered_rows
import numpy as np

def wronskian(coeffsP, coeffsQ):
    import sympy as sp
    t=sp.symbols('t')
    P=sum(c*t**i for i,c in enumerate(coeffsP))
    Q=sum(c*t**i for i,c in enumerate(coeffsQ))
    W=sp.expand(P*sp.diff(Q,t)-sp.diff(P,t)*Q)
    poly=sp.Poly(W,t)
    deg=poly.degree()
    if deg<0:
        return []
    return [int(poly.nth(i)) for i in range(deg+1)], W

for j in [2,3]:
    rows={}
    for n,c in ordered_rows(j,10):
        rows[n]=c
    for n in [3,4,5,6]:
        res=wronskian(rows[n],rows[n+1])
        if not res: 
            print(j,n,"zero"); continue
        W,_=res
        # coeffs ascending -> descending for np.roots
        desc=W[::-1]
        rts=np.roots(desc)
        print(f"j={j} n={n} deg={len(W)-1}")
        for r in rts:
            print(f"   {r.real:.6f} {r.imag:+.6f}i")
