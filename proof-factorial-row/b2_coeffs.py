"""Coefficients of -B2."""
import sympy as sp
from explore import T
xs = sp.symbols('x')
for r in range(5, 12):
    A = T[r-1].as_expr(); C = T[r-2].as_expr()
    Ap, Cp = sp.diff(A, xs), sp.diff(C, xs)
    B2 = xs**2*(3*A*Cp - 2*Ap*C) + xs**3*(A*sp.diff(Cp, xs) - Cp*Ap)
    P = sp.Poly(sp.expand(-B2), xs)
    print('r=%d -B2 coeffs asc:' % r, [int(v) for v in reversed(P.all_coeffs())])
