"""Check W^{(j)}_r = -W(U^j F_{r-1}, x U^j F_{r-2}) coefficient positivity for j=0..3."""
import sympy as sp
from explore import T
xs = sp.symbols('x')
Uop = lambda f: sp.expand(2*f + xs*sp.diff(f, xs))
for r in range(5, 11):
    out = []
    for j in range(0, 4):
        p = Uop(T[r-1].as_expr()) if False else T[r-1].as_expr()
        q = T[r-2].as_expr()
        for _ in range(j):
            p, q = Uop(p), Uop(q)
        R = xs*q
        W = sp.expand(sp.diff(p,xs)*R - p*sp.diff(R,xs))
        P = sp.Poly(-W, xs)
        coeffs = [int(v) for v in reversed(P.all_coeffs())]
        out.append(all(c >= 0 for c in coeffs))
    print('r=%d positive coefficients of -W^(j) for j=0..3:' % r, out)
