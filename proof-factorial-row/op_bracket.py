"""Check B/delta sign for p=F_{r-2}, R=x F_{r-3}."""
import sympy as sp
from explore import T
xs = sp.symbols('x')
for r in range(6, 11):
    p = T[r-2].as_expr(); q = T[r-3].as_expr(); R = xs*q
    a = sp.diff(p,xs)/p; u = sp.diff(R,xs)/R
    d = a-u
    up = sp.diff(u,xs); dp = sp.diff(d,xs)
    B = d**2*xs*(u*xs+1) + d*(u**2*xs**2+3*u*xs-up*xs**2+3) + dp*xs*(u*xs+1) - (2*u*up*xs**2+u+3*up*xs)
    ratio = sp.simplify(B/d)
    # sample points
    vals = [float(ratio.subs(xs,v)) for v in [-100,-10,-3,-1,-0.3,-0.1,-0.01,0.01,1]]
    print('r=%d B/delta signs:' % r, [1 if v>0 else -1 for v in vals])
