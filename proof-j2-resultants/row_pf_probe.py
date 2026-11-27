"""Probe: for fixed n,x, is (a_{n,r}(x))_r a PF sequence (Toeplitz matrix TN)?
Equivalent to O_n(x,t) having all real nonpositive roots.
"""
from __future__ import annotations

import itertools
import sympy as sp

x, t = sp.symbols("x t")


def rows_x(N):
    out = [sp.Integer(1)]
    for n in range(1, N + 1):
        p = n * x * out[n - 1]
        p += t * sum(sp.binomial(n, r) * out[n - r] for r in range(2, n + 1))
        out.append(sp.expand(p))
    return out


def check_pf(vals):
    D = len(vals) - 1
    N = D + 1
    T = [[sp.Integer(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            if i >= j:
                T[i][j] = vals[i - j]  # lower triangular Toeplitz
    worst = sp.Integer(0)
    for size in range(2, N + 1):
        for rows in itertools.combinations(range(N), size):
            for cols in itertools.combinations(range(N), size):
                M = sp.Matrix([[T[i][j] for j in cols] for i in rows])
                det = M.det()
                if det < worst:
                    worst = det
    return worst


def main():
    N = 9
    O = rows_x(N + 1)
    for n in range(2, N + 1):
        D = n // 2
        vals_sym = [sp.expand(O[n].coeff(t, r)) for r in range(D + 1)]
        for xv in [sp.Rational(1, 4), sp.Integer(1), sp.Integer(4)]:
            vals = [sp.N(v.subs(x, xv), 40) for v in vals_sym]
            worst = check_pf(vals)
            print(f"n={n} x={xv}: min Toeplitz minor = {float(worst):.4g}", flush=True)


if __name__ == "__main__":
    main()
