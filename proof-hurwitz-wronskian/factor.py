import sys, math
sys.path.insert(0, 'C:/Users/Abdelli/Documents/ordered-big-blocks')
from scan import ordered_rows
import sympy as sp

def wronskian(coeffsP, coeffsQ):
    t=sp.symbols('t')
    P=sum(c*t**i for i,c in enumerate(coeffsP))
    Q=sum(c*t**i for i,c in enumerate(coeffsQ))
    W=sp.expand(P*sp.diff(Q,t)-sp.diff(P,t)*Q)
    return W

for j in [2,3]:
    rows={}
    for n,c in ordered_rows(j,12):
        rows[n]=c
    for n in [5,6]:
        W=wronskian(rows[n],rows[n+1])
        print(f"j={j} n={n}")
        print("  factor W:", sp.factor(W))
        t=sp.symbols('t')
        # even/odd
        pe = sp.Poly(W, t).nth(0)  # no
        We = sum(sp.Poly(W,t).nth(k) for k in range(0, sp.Poly(W,t).degree()+1,2))*0
        # build even polynomial in u=t^2
        u=sp.symbols('u')
        pW=sp.Poly(W,t)
        deg=pW.degree()
        even=sum(pW.nth(k)*u**(k//2) for k in range(0,deg+1,2))
        odd=sum(pW.nth(k)*u**((k-1)//2) for k in range(1,deg+1,2))
        print("  even(u)=", sp.factor(even))
        print("  odd(u) =", sp.factor(odd))
PY