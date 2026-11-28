"""Test generic two-row model with c parameter (evidence / possible obstruction)."""
import sympy as sp
t,c=sp.symbols('t c')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
for u in [(1,3),(1,2,5),(1,3,7)]:
    S=1
    for x in u: S*= (1+x*t)
    A=E(S); B=t*U(S)
    P=sp.expand(A+c*B); T=sp.expand(E(sp.diff(P,t)))
    H=sp.Poly(sp.expand(W(P,T)),t)
    Hc=H.as_expr()
    print('u=',u,'degP',sp.degree(P,t))
    for cv in [0,1,2,3,5,10,50]:
        hc=sp.Poly(Hc.subs(c,cv),t)
        print('   c=',cv,'real roots',hc.count_roots(-sp.oo,sp.oo),'coeff sign min',min(hc.all_coeffs()))
