"""Exact obstruction: nonnegative Riccati coefficients are not sufficient.

For G'=t^2 z G^2+G, the coefficient O_2 is 1+t^2, which has no real root.
"""
from __future__ import annotations

from math import comb
import json
import sympy as sp
from exact_algebra import add, mul, scale


def H_row(rows, n):
    out = [0]
    for r in range(n + 1):
        out = add(out, scale(mul(rows[r], rows[n - r]), comb(n, r)))
    return out


def riccati_rows(q, nmax):
    rows = [[1]]
    hs = []
    for n in range(nmax):
        hs.append(H_row(rows, n))
        nxt = rows[n]
        for r in range(1, n + 1):
            if r < len(q) and q[r] != [0]:
                nxt = add(nxt, scale(mul(q[r], hs[n - r]), comb(n, r)))
        rows.append(nxt)
    return rows


def main():
    # q_r is the coefficient of z^r/r! in the potential Q.
    q = [[0], [0, 0, 1]]  # Q=t^2 z
    rows = riccati_rows(q, 2)
    t = sp.Symbol("t")
    p = sp.Poly(sum(sp.Integer(c) * t**i for i, c in enumerate(rows[2])), t)
    print(json.dumps({
        "potential": "t^2 z",
        "O_0": rows[0],
        "O_1": rows[1],
        "O_2": rows[2],
        "discriminant": str(p.discriminant()),
        "real_root_count": int(p.count_roots(-sp.oo, sp.oo)),
    }, indent=2))


if __name__ == "__main__":
    main()
