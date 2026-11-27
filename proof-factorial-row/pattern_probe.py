"""Detailed root patterns of F_{r-1}, U_r, xU_r, F_r."""
import numpy as np
import sympy as sp
from explore import T, U
xs = sp.symbols('x')


def rr(p):
    if p.degree() == 0:
        return np.array([])
    r = np.roots([float(v) for v in p.all_coeffs()])
    assert max(abs(r.imag)) < 1e-5, p
    return np.sort(r.real)


for r in range(5, 11):
    F = T[r]
    A = T[r - 1]
    B = T[r - 2]
    U_r = U[r]
    P = sp.Poly(sp.expand(xs * U_r.as_expr()), xs)
    pattern = sorted([(v, 'A') for v in rr(A)] +
                     [(v, 'u') for v in rr(U_r)] +
                     [(v, 'P') for v in rr(P)] +
                     [(v, 'F') for v in rr(F)])
    print('r=%d' % r, ''.join(s for _, s in pattern))
    # also raw values for small r
    if r <= 8:
        print('   A:', rr(A), ' U:', rr(U_r), ' P:', rr(P), ' F:', rr(F))
