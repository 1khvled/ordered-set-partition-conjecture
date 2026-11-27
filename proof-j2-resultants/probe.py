"""Structural probe: sign pattern of Disc/Res as polynomials in x (or s=1/x).

Exploratory only; not a proof.
"""
from __future__ import annotations

import sympy as sp

x, t, s = sp.symbols("x t s")


def rows_x(N: int):
    out = [sp.Integer(1)]
    for n in range(1, N + 1):
        p = n * x * out[n - 1]
        p += t * sum(sp.binomial(n, r) * out[n - r] for r in range(2, n + 1))
        out.append(sp.expand(p))
    return out


def sign_pattern(poly, var):
    p = sp.Poly(sp.expand(poly), var)
    c = p.all_coeffs()  # descending
    k = 0
    while k < len(c) and c[k] == 0:
        k += 1
    lead = c[k:]
    j = len(lead) - 1
    while j >= 0 and lead[j] == 0:
        j -= 1
    core = lead[: j + 1]  # strip trailing zeros
    signs = set(sp.sign(v) for v in core)
    return signs, len(core), core


def main():
    N = 9
    O = rows_x(N + 1)
    for n in range(2, N + 1):
        d = sp.expand(sp.discriminant(O[n], t))
        signs, ln, core = sign_pattern(d, x)
        print(f"Disc n={n}: deg_x={sp.Poly(d, x).degree()} signs={signs} "
              f"len={ln} c0={core[0]} c_last={core[-1]}")

    for n in range(1, N):
        r = sp.expand(sp.resultant(O[n], O[n + 1], t))
        signs, ln, core = sign_pattern(r, x)
        print(f"Res n={n},{n+1}: deg_x={sp.Poly(r, x).degree()} signs={signs} "
              f"len={ln} c0={core[0]} c_last={core[-1]}")

    # Also look in s=1/x variable for Disc_t(Q) where Q_n(s,t)=O_n(1/s,t)*s^n
    print("--- s variable ---")
    for n in range(2, N + 1):
        Q = sp.expand(s**n * O[n].subs(x, 1 / s))
        d = sp.expand(sp.discriminant(Q, t))
        signs, ln, core = sign_pattern(d, s)
        print(f"sDisc n={n}: deg_s={sp.Poly(d, s).degree()} signs={signs} "
              f"len={ln} c0={core[0]} c_last={core[-1]}")


if __name__ == "__main__":
    main()
