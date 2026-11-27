"""Fit -W_r as positive combination of F-squares and previous W's with x^{2j}."""
import sympy as sp
from explore import T
xs = sp.symbols('x')


def W(r):
    F = T[r].as_expr()
    G = xs * T[r - 1].as_expr()
    return sp.expand(sp.diff(F, xs) * G - F * sp.diff(G, xs))


for r in range(5, 13):
    target = sp.Poly(-W(r), xs)
    names, polys = [], []
    for j in range(0, 4):
        for t in (r-1, r-2, r-3):
            if t >= 2:
                names.append('x%d F%d^2' % (2*j, t))
                polys.append(xs**(2*j)*T[t].as_expr()**2)
        for t in (r-1, r-2, r-3):
            if t >= 4:
                names.append('x%d (-W%d)' % (2*j, t))
                polys.append(xs**(2*j)*(-W(t)))
    N = target.degree()+1
    M = sp.zeros(N, len(polys)); rhs = sp.zeros(N, 1)
    for d in range(N):
        rhs[d] = target.coeff_monomial(xs**d)
        for i, p in enumerate(polys):
            M[d, i] = sp.Poly(sp.expand(p), xs).coeff_monomial(xs**d)
    sol = sp.linsolve((M, rhs))
    if not sol:
        print('r=%d no exact solution' % r); continue
    ss = list(sol)[0]
    nz = [(nm, v) for nm, v in zip(names, ss) if v != 0]
    print('r=%d:' % r, nz[:12])
