"""Inspect W(T_n, Q) and its recurrence pieces."""
import sympy as sp
from common import Family, t
F=Family(20); E,W=Family.E,Family.W
for n in (8,10,12):
    Pn=F.Qpoly(n); Q=F.Qpoly(n-2)+E(F.Qpoly(n-2))  # U P_{n-2}
    Tn=F.T(Pn); Tn1=F.T(F.Qpoly(n-1))
    w=W(Tn,Q)
    p1=W(E(Tn1)+Tn1,Q); p2=W(E(E(Q)),Q)
    print('n=',n,'W(Tn,Q) roots',sp.Poly(w,t).count_roots(-sp.oo,sp.oo),'lc',sp.LC(w,t))
    print('   W((1+E)T_{n-1},Q) roots',sp.Poly(p1,t).count_roots(-sp.oo,sp.oo),'lc',sp.LC(p1,t))
    print('   W(E^2Q,Q) roots',sp.Poly(p2,t).count_roots(-sp.oo,sp.oo),'lc',sp.LC(p2,t))
# generic Q test
Q=1+6*t+11*t**2+6*t**3
print('generic W(E^2Q,Q)=',sp.factor(W(E(E(sp.Poly(Q,t))),sp.Poly(Q,t)).as_expr()))
