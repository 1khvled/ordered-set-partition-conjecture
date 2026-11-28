"""Wronskian N = W(d_x psi_n, psi_n') : coefficient signs and factorization."""
import sympy as sp
from probe2 import psi

x, u = sp.symbols('x u')
for n in range(3, 11):
    f = psi(n).as_expr()
    dPdx = sp.diff(f, x)
    dPdu = sp.diff(f, u)
    N = sp.expand(sp.diff(dPdx, u)*dPdu - dPdx*sp.diff(f, u, 2))
    coeffs = sp.Poly(N, u).all_coeffs()
    N1 = sp.expand(N.subs(x, 1))
    coeffs1 = sp.Poly(N1, u).all_coeffs()
    allneg = all(c <= 0 for c in coeffs1)
    print('n=%d degN=%d all-neg coeffs at x=1? %s' % (n, sp.Poly(N,u).degree(), allneg))
    print('   N(x=1)=', sp.factor(N1) if n<=8 else str(N1)[:250]+'...')
