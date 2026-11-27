import sys, math
sys.path.insert(0, 'C:/Users/Abdelli/Documents/ordered-big-blocks')
from scan import ordered_rows
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

j=2
rows={}
for n,c in ordered_rows(j,16):
    rows[n]=c
for n in range(0,15):
    W=wronskian(rows[n],rows[n+1])
    print(f"n={n} deg={len(W)-1 if W else None} W={W}")
PY