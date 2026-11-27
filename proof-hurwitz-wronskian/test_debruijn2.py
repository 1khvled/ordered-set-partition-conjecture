import numpy as np
import sympy as sp

def W_roots(P,Q):
    t=sp.symbols('t')
    p=1; q=1
    for r in P: p*=(t-r)
    for r in Q: q*=(t-r)
    W=sp.expand(p*sp.diff(q,t)-sp.diff(p,t)*q)
    poly=sp.Poly(W,t); deg=poly.degree()
    if deg<0: return []
    cs=[int(poly.nth(i)) for i in range(deg+1)]
    return np.roots(cs[::-1])

# same degree interlacing: P roots -4,-2 ; Q roots -3,-1
P=[-4,-2]; Q=[-3,-1]
print(W_roots(P,Q))
# another same degree
P=[-5,-3,-1]; Q=[-4,-2,0]  # last root 0 not negative but test
print(W_roots(P,Q))
PY