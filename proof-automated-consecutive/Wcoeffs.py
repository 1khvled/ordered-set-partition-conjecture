import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs
import sympy as sp

t = sp.Symbol('t')

def Wcoeffs(j,nmax):
    O={}
    for n,coeffs in ordered_rows(j,nmax):
        O[n]=poly_from_coeffs(coeffs)
    for n in range(1,nmax):
        P=O[n]; Q=O[n+1]
        if P.degree()<=0: continue
        W=sp.expand(P.as_expr()*sp.diff(Q.as_expr(),t) - sp.diff(P.as_expr(),t)*Q.as_expr())
        Wp=sp.Poly(W,t)
        coeffs=Wp.all_coeffs()[::-1]  # ascending
        neg=[c for c in coeffs if c<0]
        print(f"j={j} n={n} deg={Wp.degree()} #negcoeffs={len(neg)} mincoeff={min(coeffs)}")

for j in [2,3,4]:
    Wcoeffs(j,12)