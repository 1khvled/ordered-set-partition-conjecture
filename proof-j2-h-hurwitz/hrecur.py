"""Look for a simple recurrence for H_n itself."""
import sympy as sp
from common import Family, t
F=Family(20)
for n in range(8,16):
    Hn=F.H(n); Hn1=F.H(n-1); Hn2=F.H(n-2)
    print('n',n,'deg',Hn.degree(),'degprev',Hn1.degree())
    combos={'H_n - t^2 H_{n-2}':Hn - sp.Poly(t,t)**2*Hn2,
            'H_n - t H_{n-1}*?':None}
    for name,c in list(combos.items()):
        if c is None: continue
        print('  ',name,'realroots',sp.Poly(c.as_expr(),t).count_roots(-sp.oo,sp.oo),'coeffmin',min(c.all_coeffs()))
