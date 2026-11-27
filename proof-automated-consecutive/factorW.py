import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs
import sympy as sp

t = sp.Symbol('t')

def factor_W(j,nmax):
    O={}
    for n,coeffs in ordered_rows(j,nmax):
        O[n]=poly_from_coeffs(coeffs)
    for n in range(1,nmax):
        P=O[n]; Q=O[n+1]
        if P.degree()<=0: continue
        W=sp.expand(P.as_expr()*sp.diff(Q.as_expr(),t) - sp.diff(P.as_expr(),t)*Q.as_expr())
        Wp=sp.Poly(W,t)
        print(f"j={j} n={n} degW={Wp.degree()} factor={sp.factorint(Wp.as_expr()) if False else 'NA'}")
        print("   W=",Wp.as_expr())
        # try factor over QQ
        f=sp.factor(Wp.as_expr())
        print("   factored:", f)

factor_W(2,6)