"""Euclidean/continued-fraction data for T_n/P_n."""
import sympy as sp
from common import Family, t
F=Family(20)
for n in (6,7,8,9,10,11,12):
    P=F.Qpoly(n).as_expr(); Tf=F.T(F.Qpoly(n)).as_expr()
    q,r=sp.div(sp.Poly(P,t),sp.Poly(Tf,t))
    print('n=',n,'degP',sp.degree(P,t),'degT',sp.degree(Tf,t))
    print('  quotient',sp.factor(q.as_expr()),'remainder deg',sp.degree(r.as_expr(),t),'lc',sp.LC(r,t))
    # next step
    q2,r2=sp.div(Tf,r.as_expr(),t)
    print('  q2',sp.factor(q2),'r2',sp.factor(r2))
