"""Search for an inductive positivity certificate for -W_r as a nonneg combination
of inductively-positive building blocks (squares and -W_{r-1})."""
import sympy as sp
from explore import T
xs = sp.symbols('x')


def W(r):
    F = T[r].as_expr()
    G = xs * T[r - 1].as_expr()
    return sp.expand(sp.diff(F, xs) * G - F * sp.diff(G, xs))


for r in range(5, 13):
    target = sp.Poly(-W(r), xs)
    basis = []
    names = []
    for j in range(0, 4):
        basis.append(sp.Poly(sp.expand(xs ** (2*j) * (-W(r-1))), xs))
        names.append('x^{2j}*(-W_{r-1})' % () if False else 'x^%d*(-W_%d)' % (2*j, r-1))
    for s in [r-1, r-2, r-3]:
        if s < 2:
            continue
        f = T[s].as_expr()
        for j in range(0, 5):
            basis.append(sp.Poly(sp.expand(xs ** (2*j) * f ** 2), xs))
            names.append('x^%d*F_%d^2' % (2*j, s))
    # build coefficient matrix on monomials up to degree of target
    N = target.degree() + 1
    rows, rhs = [], []
    for d in range(N):
        rows.append([int(b.coeff_monomial(xs ** d)) for b in basis])
        rhs.append(int(target.coeff_monomial(xs ** d)))
    M = sp.Matrix(rows); v = sp.Matrix(rhs)
    sol = sp.linsolve((M, v))
    if not sol:
        print('r=%d no solution' % r); continue
    s = list(sol)[0]
    print('r=%d solution (free params): %s' % (r, s))
