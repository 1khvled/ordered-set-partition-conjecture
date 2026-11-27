import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs
import sympy as sp

t = sp.Symbol('t')

def Wpos(j,nmax):
    O={}
    for n,coeffs in ordered_rows(j,nmax):
        O[n]=poly_from_coeffs(coeffs)
    for n in range(1,nmax):
        P=O[n]; Q=O[n+1]
        if P.degree()<=0: continue
        W=sp.expand(P.as_expr()*sp.diff(Q.as_expr(),t) - sp.diff(P.as_expr(),t)*Q.as_expr())
        Wp=sp.Poly(W,t)
        # check real roots
        roots=sp.nroots(Wp.as_expr()) if Wp.degree()>0 else []
        realroots=[r for r in roots if abs(sp.im(r))<1e-30]
        lc=sp.Poly(W,t).nth(Wp.degree())
        print(f"j={j} n={n} degW={Wp.degree()} #realroots={len(realroots)} lc_sign={sp.sign(lc)}")
        if len(realroots)>0:
            print("   real roots:", [complex(r) for r in realroots[:5]])

for j in [2,3,4]:
    Wpos(j,12)