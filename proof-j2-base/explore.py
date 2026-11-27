"""Small exact/numerical probes for the j=2 base case.

No assertion in this file is an all-index proof.  It searches for a sharper
rootwise quantity than the consecutive residue itself.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb


def rows(N: int) -> list[list[int]]:
    """Ascending integer coefficients of O_{n,2}, via the first-block rule."""
    out = [[1]]
    for n in range(1, N + 1):
        prev = out[n - 1]
        row = [n * x for x in prev] + [0]
        for r in range(2, n + 1):
            old = out[n - r]
            c = comb(n, r)
            for k, x in enumerate(old):
                row[k + 1] += c * x
        while row and row[-1] == 0:
            row.pop()
        out.append(row)
    return out


def peval(a: list[int], x):
    y = 0
    for c in reversed(a):
        y = y * x + c
    return y


def pder(a: list[int]) -> list[int]:
    return [k * a[k] for k in range(1, len(a))]


def convolution(a: list[int], b: list[int]) -> list[int]:
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def main() -> None:
    import sympy as sp

    x = sp.symbols("x")
    R = rows(18)
    worst = None
    for n in range(2, len(R)):
        p = sum(c * x**k for k, c in enumerate(R[n]))
        u = [0] * n
        for k in range(n):
            prod = convolution(R[k], R[n - 1 - k])
            for j, c in enumerate(prod):
                u[j] += comb(n - 1, k) * c
        q = sum(c * x**j for j, c in enumerate(u))
        roots = sp.nroots(p, n=50, maxsteps=200)
        vals = []
        for mu in roots:
            if abs(sp.im(mu)) > sp.Float("1e-40"):
                vals.append(("nonreal", mu))
                continue
            mu = sp.re(mu)
            ratio = sp.N(q.subs(x, mu) / sp.diff(p, x).subs(x, mu), 40)
            target = sp.N(2 * mu + n * (1 + mu) * ratio, 40)
            vals.append((sp.N(mu, 20), ratio, target))
            if worst is None or ratio < worst[0]:
                worst = (ratio, n, mu, target)
        print(f"n={n} deg={len(R[n])-1} U/O'={vals}")
    print("worst", worst)


if __name__ == "__main__":
    main()
