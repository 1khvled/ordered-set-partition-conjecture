"""Compute scaled derivative polynomials D_n = x d_x psi_n and look for closed form."""
import sympy as sp
from probe2 import psi

x, u = sp.symbols('x u')
for n in range(2, 9):
    Ps = psi(n)
    D = sp.expand(x * sp.diff(Ps.as_expr(), x))
    print('n=%d psi=%s' % (n, sp.factor(Ps.as_expr())))
    print('    D=%s' % sp.factor(D))
