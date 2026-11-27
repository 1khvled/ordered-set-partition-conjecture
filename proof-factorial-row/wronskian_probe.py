"""Exact Wronskians W(F_r,G_r) and factorizations."""
import sympy as sp
from explore import T, U, MAXR
xs = sp.symbols('x')

for r in range(4, 11):
    F = T[r].as_expr()
    G = xs * T[r - 1].as_expr()
    W = sp.factor(sp.expand(sp.diff(F, xs) * G - F * sp.diff(G, xs)))
    Fprev = T[r - 1].as_expr()
    Gprev = xs * T[r - 2].as_expr()
    Wp = sp.factor(sp.expand(sp.diff(Fprev, xs) * Gprev - Fprev * sp.diff(Gprev, xs)))
    print('r=%d' % r)
    print('  W_r =', W)
    print('  W_r/W_{r-1} =', sp.factor(W / Wp) if Wp != 0 else None)
    print('  W_r / (x^? F?) ...', sp.factor(W / (xs ** (r - 4) * Fprev)))
