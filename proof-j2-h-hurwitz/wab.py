"""Inspect W(A,B) for the two-term recurrence A=E P_{n-1}, B=t P_{n-2,a+1}."""
import sympy as sp
from common import Family, t
F=Family(20); E,T,W=Family.E,Family.T,Family.W
for a in (1,2):
    print('=== a',a)
    for n in (8,10,12):
        A=E(F.Qpoly(n-1,a)); B=sp.Poly(t,t)*F.Qpoly(n-2,a+1)
        w=W(A,B)
        print(' n',n,'W(A,B) real roots',sp.Poly(w,t).count_roots(-sp.oo,sp.oo),'lc',sp.LC(w,t),'coeffmin',min(w.all_coeffs()))
        print('   factor',sp.factor(w.as_expr()) if n<=10 else '')
