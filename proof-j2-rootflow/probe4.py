"""Compare roots of A_n = -d_x psi_n/u with roots of psi_n."""
import sympy as sp
import mpmath as mp
from probe2 import psi

x, u = sp.symbols('x u')


def A_poly(n, X):
    Ps = psi(n)
    D = sp.diff(Ps.as_expr(), x)
    Xu = sp.Rational(X)
    return sp.Poly(sp.cancel(-D.subs(x, Xu) / u), u)


def roots_of(P):
    cs = [mp.mpf(str(c)) for c in P.all_coeffs()]
    return sorted([mp.re(r) for r in mp.polyroots(cs, maxsteps=300, error=False)])


mp.mp.dps = 40
for n in range(4, 11):
    P = psi(n, sp.Rational(1))
    A = A_poly(n, 1)
    rp, ra = roots_of(P), roots_of(A)
    print('n=%d deg psi=%d deg A=%d' % (n, P.degree(), A.degree()))
    print('  psi roots:', [round(float(v), 5) for v in rp])
    print('  A   roots:', [round(float(v), 5) for v in ra])
