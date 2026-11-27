import sys, math
sys.path.insert(0, 'C:/Users/Abdelli/Documents/ordered-big-blocks')
from scan import ordered_rows
import numpy as np
import sympy as sp

def wronskian(coeffsP, coeffsQ):
    t=sp.symbols('t')
    P=sum(c*t**i for i,c in enumerate(coeffsP))
    Q=sum(c*t**i for i,c in enumerate(coeffsQ))
    W=sp.expand(P*sp.diff(Q,t)-sp.diff(P,t)*Q)
    poly=sp.Poly(W,t)
    deg=poly.degree()
    if deg<0: return None
    return [int(poly.nth(i)) for i in range(deg+1)]

for j in [2,3]:
    rows={}
    for n,c in ordered_rows(j,12):
        rows[n]=c
    for n in [5,6,7]:
        W=wronskian(rows[n],rows[n+1])
        if not W: continue
        deg=len(W)-1
        # even part coeffs for t^2: W_even[k] = W[2k]; odd: W[2k+1]
        even=[W[2*k] for k in range((deg//2)+1)]
        odd=[W[2*k+1] for k in range((deg//2))]  # t * sum odd[k] t^{2k}
        print(f"j={j} n={n} deg={deg}")
        for name,cs in [("even",even),("odd",odd)]:
            desc=cs[::-1]
            rts=np.roots(desc) if len(desc)>1 else []
            print(f"  {name} deg={len(cs)-1} roots:")
            for r in rts:
                print(f"     {r.real:.6f} {r.imag:+.6f}i")
PY