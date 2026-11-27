import sys, math
sys.path.insert(0, 'C:/Users/Abdelli/Documents/ordered-big-blocks')
from scan import ordered_rows

def wronskian(coeffsP, coeffsQ):
    import sympy as sp
    t=sp.symbols('t')
    P=sum(c*t**i for i,c in enumerate(coeffsP))
    Q=sum(c*t**i for i,c in enumerate(coeffsQ))
    W=sp.expand(P*sp.diff(Q,t)-sp.diff(P,t)*Q)
    poly=sp.Poly(W,t)
    deg=poly.degree()
    if deg<0:
        return [], W
    cs=[int(poly.nth(i)) for i in range(deg+1)]
    return cs, W

for j in [2,3,4]:
    rows={}
    for n,c in ordered_rows(j,12):
        rows[n]=c
    for n in range(0,8):
        W,expr=wronskian(rows[n],rows[n+1])
        print(f"j={j} n={n} degP={len(rows[n])-1} degQ={len(rows[n+1])-1} degW={len(W)-1 if W else None}")
        print("  W coeffs asc:", W)
PY