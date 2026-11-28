"""Signs of recurrence pieces at roots of Q_n."""
import numpy as np, sympy as sp
from common import Family, t
F=Family(20); E=Family.E
for n in (10,12):
    Pn=F.Qpoly(n); Q=F.Qpoly(n-2)+E(F.Qpoly(n-2))
    Tn=F.T(Pn); Tn1=F.T(F.Qpoly(n-1)); c=n-1
    part1=sp.expand(Tn1.as_expr()+t*sp.diff(Tn1.as_expr(),t)+Tn1.as_expr())  # U T_{n-1}
    part2=sp.expand(E(E(Q.as_expr())).as_expr())*c
    qroots=[r for r in sp.nroots(Q, n=20) if abs(sp.im(r))<sp.Rational(1,10**12)]
    print('n',n)
    for qr in qroots:
        x=float(sp.re(qr))
        v1=float(part1.subs(t,x)); v2=float(part2.subs(t,x))
        print('  mu=',x,' U T_{n-1}=',v1,' c E^2Q=',v2,' sum=',v1+v2)
