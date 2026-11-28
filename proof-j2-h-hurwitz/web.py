"""Root web for A,B,alpha,beta."""
import sympy as sp, numpy as np
from common import Family, t
F=Family(20); E=Family.E
def roots(p):
    cs=[float(c) for c in p.all_coeffs()]
    return sorted([z.real for z in np.roots(cs) if abs(z.imag)<1e-8])
n=10; a=1
S=F.Qpoly(n-1,a); K=F.Qpoly(n-2,a+1); Tprev=F.T(S)
A=E(S); B=sp.Poly(t,t)*K
alpha=sp.Poly(sp.expand(Tprev.as_expr()+E(Tprev.as_expr()).as_expr()),t)
beta=E(E(K).as_expr())
for name,p in [('A=E S',A),('B=tK',B),('alpha=U Tprev',alpha),('beta=E^2K',beta)]:
    print(name, ['%.6g'%x for x in roots(p)])
