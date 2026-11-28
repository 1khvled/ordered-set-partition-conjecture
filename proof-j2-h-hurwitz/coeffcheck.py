"""Coefficient inequalities for H_n."""
import sympy as sp
from common import Family
F=Family(16)
for n in (8,10,12,14):
    c=list(reversed(F.H(n).all_coeffs()))  # ascending
    d=len(c)-1
    lc=[sp.simplify(c[k]**2-c[k-1]*c[k+1])>0 for k in range(1,d)]
    newt=[sp.simplify(c[k]**2 - sp.Rational(k+1,k)*sp.Rational(d-k+1,d-k)*c[k-1]*c[k+1])>0 for k in range(1,d)]
    print('n',n,'logconc',all(lc),'newton',all(newt), 'ratios', [float(c[k+1]/c[k]) for k in range(d)])
