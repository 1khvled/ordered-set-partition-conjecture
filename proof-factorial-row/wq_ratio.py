"""Compare W(D_p, x D_q) with W(p, x q)."""
import sympy as sp
from explore import T
xs = sp.symbols('x')
for r in range(6, 13):
    p = T[r-2].as_expr(); q = T[r-3].as_expr()
    Pp = 2*p + xs*sp.diff(p, xs)
    Q = xs*(2*q + xs*sp.diff(q, xs))
    Wnew = sp.expand(sp.diff(Pp,xs)*Q - Pp*sp.diff(Q,xs))
    Wold = sp.expand(sp.diff(p,xs)*(xs*q) - p*sp.diff(xs*q,xs))
    ratio = sp.factor(sp.cancel(Wnew/Wold))
    print('r=%d  W_old=%s' % (r, sp.factor(Wold)))
    print('   W_new/W_old =', ratio)
