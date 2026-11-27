"""Evaluate B2 at roots of A=F_{r-1} and C=F_{r-2}, and globally via Sturm/sampling."""
import sympy as sp
import numpy as np
from explore import T
xs = sp.symbols('x')

for r in range(5, 11):
    A = T[r-1].as_expr(); C = T[r-2].as_expr()
    Ap = sp.diff(A, xs); App = sp.diff(Ap, xs)
    Cp = sp.diff(C, xs); Cpp = sp.diff(Cp, xs)
    B2 = xs**2*(3*A*Cp - 2*Ap*C) + xs**3*(A*Cpp - Cp*Ap)
    ra = [float(v) for v in np.sort(np.roots([float(v) for v in sp.Poly(A, xs).all_coeffs()]))]
    rc = [float(v) for v in np.sort(np.roots([float(v) for v in sp.Poly(C, xs).all_coeffs()]))]
    va = [float(sp.N(B2.subs(xs, sp.Float(v, 30)))) for v in ra]
    vc = [float(sp.N(B2.subs(xs, sp.Float(v, 30)))) for v in rc]
    print('r=%d A roots %s -> B2 %s' % (r, [round(v,4) for v in ra], ['%.4g'%v for v in va]))
    print('r=%d C roots %s -> B2 %s' % (r, [round(v,4) for v in rc], ['%.4g'%v for v in vc]))
