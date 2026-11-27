"""Probe: m_{n,r}(x) = deformation multiplier applied to the Fibonacci limit.
Check whether (m_{n,r}(x))_r is a PF (Pólya frequency) finite sequence for x>0:
all Toeplitz minors >= 0.  If yes, it preserves real-rootedness of Phi_n.
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


def minors_nonneg(vals):
    D = len(vals) - 1
    N = D + 1
    T = [[sp.Float(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            if j >= i:
                T[i][j] = vals[j - i]
    for size in range(1, N + 1):
        for rows in itertools.combinations(range(N), size):
            for cols in itertools.combinations(range(N), size):
                M = [[T[i][j] for j in cols] for i in rows]
                det = sp.Matrix(M).det()
                if det < -sp.Float("1e-8"):
                    return False, (size, rows, cols, float(det))
    return True, None


def main():
    N = 9
    O = rows_x(N + 1)
    for n in range(2, N + 1):
        D = n // 2
        for xv in [sp.Rational(1, 10), sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(10)]:
            vals = []
            for r in range(D + 1):
                ar = sp.expand(O[n].coeff(t, r)).subs(x, xv)
                m = 2**r * ar / (sp.factorial(n) * sp.binomial(n - r, r) * xv**(n - 2 * r))
                vals.append(sp.N(m, 30))
            ok, bad = minors_nonneg(vals)
            if not ok:
                print(f"n={n} x={xv}: PF FAILS {bad}", flush=True)
                break
        else:
            print(f"n={n}: PF holds at all tested x", flush=True)


if __name__ == "__main__":
    main()
