"""Test Bezoutian-style Cauchy-Schwarz: h1^2 <= 4 h0 h2 pointwise?"""
import sympy as sp
from common import Family, t
F=Family(16); E,W=Family.E,Family.W
for a in (1,2):
  for n in (8,10,12):
    S=F.Qpoly(n-1,a); K=F.Qpoly(n-2,a+1)
    A=E(S); B=sp.Poly(t,t)*K
    alpha=sp.Poly(sp.expand(E(E(S.as_expr())).as_expr()+E(S.as_expr()).as_expr()),t) # (1+E)T_{n-1}
    # T_{n-1}=E S' ; (1+E)T = E(ES')+ES'
    Tprev=F.T(S)
    alpha=sp.Poly(sp.expand(Tprev.as_expr()+E(Tprev.as_expr()).as_expr()),t)
    beta=E(E(K).as_expr())  # E^2 K
    h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
    D=sp.Poly(sp.expand(4*h0.as_expr()*h2.as_expr()-h1.as_expr()**2),t)
    print('a',a,'n',n,'h0roots',sp.Poly(h0.as_expr(),t).count_roots(-sp.oo,sp.oo),
          'h2roots',sp.Poly(h2.as_expr(),t).count_roots(-sp.oo,sp.oo),
          'D roots',D.count_roots(-sp.oo,sp.oo),'D lc sign',sp.sign(sp.LC(D,t)),
          'D(0)',D.eval(0))
