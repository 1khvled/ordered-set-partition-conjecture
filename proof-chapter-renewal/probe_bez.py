"""Exploration only: Bezoutian / Hurwitz structure of the boundary ratio."""
import sympy as sp
from explore import rows

t, s = sp.symbols('t s')


def bezoutian(p, q):
    """Coefficient matrix of (p(t)q(s)-p(s)q(t))/(t-s) in monomial basis."""
    n = max(sp.degree(p, t), sp.degree(q, t))
    M = sp.zeros(n, n)
    expr = sp.cancel((p * q.subs(t, s) - p.subs(t, s) * q) / (t - s))
    poly = sp.Poly(sp.expand(expr), t, s)
    for (i, j), c in poly.terms():
        M[i, j] = c
    return M


if __name__ == '__main__':
    for j in [2, 3]:
        O, H = rows(j, 4 * j)
        n = 2 * j
        for d in [1, j]:
            p, q = H[n], O[n - d]
            M = bezoutian(p, q)
            print('=== j=%d n=%d d=%d  size=%d' % (j, n, d, M.rows))
            print(' M =', M)
            print(' leading principal minors:', [M[:k, :k].det() for k in range(1, M.rows + 1)])
