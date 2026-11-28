"""Factor D_n = H_n - t^2 H_{n-2}."""
import sympy as sp
from common import Family, t
F=Family(16)
for n in (6,7,8,9,10):
    Hn=F.H(n); Hn2=F.H(n-2) if n>=6 else None
    D=Hn - sp.Poly(t,t)**2*Hn2
    print('n',n,'factor=',sp.factor(D.as_expr()))
