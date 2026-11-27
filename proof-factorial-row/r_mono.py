"""Verify -B2 = x^2 (2C+xC')^2 R'(x), R=A/(2C+xC'); check R monotone via interlacing."""
import sympy as sp
import numpy as np
from explore import T
xs = sp.symbols('x')

for r in range(5, 12):
    A = T[r-1].as_expr(); C = T[r-2].as_expr()
    Ap, Cp = sp.diff(A, xs), sp.diff(C, xs)
    B2 = xs**2*(3*A*Cp - 2*Ap*C) + xs**3*(A*sp.diff(Cp, xs) - Cp*Ap)
    D = 2*C + xs*Cp
    Rp = sp.simplify(sp.diff(A/D, xs))
    ident = sp.simplify(sp.together(B2 + xs**2*D**2*Rp))
    assert ident == 0, (r, ident)
    # roots of A and D
    ra = [float(v) for v in np.sort(np.roots([float(v) for v in sp.Poly(A, xs).all_coeffs()]))]
    rd = [float(v) for v in np.sort(np.roots([float(v) for v in sp.Poly(D, xs).all_coeffs()]))]
    allr = sorted([(v,'A') for v in ra] + [(v,'D') for v in rd])
    pat = ''.join(s for _,s in allr)
    WA = sp.expand(Ap*D - A*sp.diff(D, xs))
    pts = [-100, -3, -1, -0.1, -0.01, 0.01, 1]
    signs = [1 if float(WA.subs(xs,v))>0 else -1 for v in pts]
    print('r=%d pattern %s monotone-sign(N=A\'D-AD\') %s' % (r, pat, signs))
