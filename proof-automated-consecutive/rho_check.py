import sys, math, time
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs, real_roots
import sympy as sp

t=sp.Symbol('t')

def rho_check(j,nmax):
    O={}
    for n,coeffs in ordered_rows(j,nmax):
        O[n]=poly_from_coeffs(coeffs)
    worst=0
    bad=0
    for n in range(1,nmax):
        P=O[n]; Q=O[n+1]
        if P.degree()<=0: continue
        roots=real_roots(P)
        # verify count
        if len(roots)!=P.degree():
            print(f"j={j} n={n} real roots {len(roots)} != deg {P.degree()}")
            continue
        for mu in roots:
            muN=sp.N(mu,50)
            Pp=sp.diff(P.as_expr(),t).subs(t,muN)
            Qv=Q.as_expr().subs(t,muN)
            rho=sp.N(Qv/Pp,40)
            rf=float(rho)
            if rf>worst: worst=rf
            if rf>=0:
                bad+=1
                print(f"   BAD Rho>=0 j={j} n={n} rho={rf:.6e}")
    print(f"j={j} nmax={nmax} bad={bad} worst_rho={worst:.6e}")

rho_check(2,40)