"""Inspect -W_r positivity: roots, coefficient log-concavity, grouped quadratics."""
import sympy as sp
import numpy as np
from explore import T
xs = sp.symbols('x')

for r in range(4, 13):
    F = T[r].as_expr()
    G = xs * T[r - 1].as_expr()
    W = sp.expand(sp.diff(F, xs) * G - F * sp.diff(G, xs))
    P = sp.Poly(-W, xs)
    c = [int(v) for v in P.all_coeffs()]  # descending
    rr = np.roots([float(v) for v in c]) if P.degree() > 0 else []
    imag = max(abs(rr.imag)) if len(rr) else 0.0
    # check discriminant of each consecutive triple in reverse (ascending a0,a1,a2,...)
    a = c[::-1]
    ds = [a[2*j-1]**2 - 4*a[2*j]*a[2*j-2] for j in range(1, (len(a))//2 + 1) if 2*j < len(a)]
    print('r=%d deg=%d max|Im(root)|=%.3g' % (r, P.degree(), imag),
          'quadratic discs (must be <0):', [float(d) for d in ds])
    print('   coeffs ascending:', a[:8])
