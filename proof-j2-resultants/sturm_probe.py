"""Exact check: count positive real roots of Disc_x and Res_x polynomials.

If either has a positive real root for some n, (D)/(R) is false (exact obstruction).
Uses Sturm sequences exactly over QQ.
"""
from __future__ import annotations

import sympy as sp

x, t = sp.symbols("x t")


def rows_x(N):
    out = [sp.Integer(1)]
    for n in range(1, N + 1):
        p = n * x * out[n - 1]
        p += t * sum(sp.binomial(n, r) * out[n - r] for r in range(2, n + 1))
        out.append(sp.expand(p))
    return out


def count_pos(poly):
    p = sp.Poly(sp.expand(poly), x)
    if p.degree() <= 0:
        return 0
    return p.count_roots(0, sp.oo)


def main():
    N = 11
    O = rows_x(N + 1)
    for n in range(2, N + 1):
        d = sp.discriminant(O[n], t)
        cp = count_pos(d)
        print(f"Disc n={n}: pos real roots in x = {cp}", flush=True)
    for n in range(1, N):
        r = sp.resultant(O[n], O[n + 1], t)
        cp = count_pos(r)
        print(f"Res n={n},{n+1}: pos real roots in x = {cp}", flush=True)


if __name__ == "__main__":
    main()
