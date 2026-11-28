"""Even/odd parts of H_n."""
import sympy as sp
from common import Family, t
F=Family(14)
s=sp.symbols('s')
for n in (8,10,12):
    H=F.H(n).as_expr()
    p=sp.Poly(H,t)
    A=sum(c*s**k for k,c in enumerate(reversed(p.all_coeffs())) if (p.degree()-k)%2==0)
    B=sum(c*s**k for k,c in enumerate(reversed(p.all_coeffs())) if (p.degree()-k)%2==1)
    # simpler: collect even/odd powers
    A=0;B=0
    for (k,),c in p.terms():
        if k%2==0: A+=c*s**(k//2)
        else: B+=c*s**((k-1)//2)
    print('n',n)
    print(' A=',sp.factor(A))
    print(' B=',sp.factor(B))
