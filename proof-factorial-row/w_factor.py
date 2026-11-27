"""Factorization of -W_r and check min imaginary part of roots."""
import sympy as sp
import numpy as np
from explore import T
xs = sp.symbols('x')

for r in range(4, 13):
    F = T[r].as_expr()
    G = xs * T[r - 1].as_expr()
    W = sp.factor(sp.expand(sp.diff(F, xs) * G - F * sp.diff(G, xs)))
    P = sp.Poly(-W, xs)
    if P.degree():
        roots = np.roots([float(v) for v in P.all_coeffs()])
        print('r=%d min|Im|=%.4g fact=%s' % (r, min(abs(roots.imag)), sp.factor(-W)))
    else:
        print('r=%d const %s' % (r, -W))
