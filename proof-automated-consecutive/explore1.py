import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs, real_roots
import sympy as sp

t = sp.Symbol('t')

def run(j,nmax):
    O={}
    for n,coeffs in ordered_rows(j,nmax):
        O[n]=poly_from_coeffs(coeffs)
    for n in range(1,nmax):
        P=O[n]; Q=O[n+1]
        if P.degree()<=0: continue
        roots=real_roots(P)
        print(f"j={j} n={n} degP={P.degree()} degQ={Q.degree()}")
        for mu in roots:
            muN=sp.N(mu,60)
            Pp=sp.diff(P.as_expr(),t).subs(t,muN)
            Qv=Q.as_expr().subs(t,muN)
            rho=sp.N(Qv/Pp,40)
            print(f"   mu={float(muN):.6f} rho={float(rho):.6e}")

for j in [2,3]:
    run(j,8)