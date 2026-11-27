"""Exact finite verification of the operator lemma W(U p, x U q) < 0.

This is finite evidence only; the all-r proof is in NOTES.md.
"""
from __future__ import annotations

import sympy as sp

x = sp.symbols("x")


def associated_table(max_r: int) -> list[list[int]]:
    A = [[0] * (max_r + 2) for _ in range(max_r + 1)]
    A[0][0] = 1
    for r in range(2, max_r + 1):
        for k in range(1, r // 2 + 1):
            A[r][k] = k * A[r - 1][k] + (r - 1) * A[r - 2][k - 1]
    return A


def E_poly(A: list[list[int]], r: int) -> sp.Expr:
    return sum(sp.factorial(k) * A[r][k] * x**k for k in range(r // 2 + 1))


def F_poly(A: list[list[int]], r: int) -> sp.Expr:
    if r < 2:
        return sp.Integer(0)
    return sp.cancel(E_poly(A, r) / x)


def U(f: sp.Expr) -> sp.Expr:
    return sp.expand(2 * f + x * sp.diff(f, x))


def W(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(f, x) * g - f * sp.diff(g, x))


def variations(seq: list[int]) -> int:
    signs = [1 if c > 0 else -1 for c in seq if c != 0]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def real_roots(poly: sp.Poly) -> int:
    seq = sp.sturm(poly)
    return variations([int(c) for c in seq[0].all_coeffs()]) - variations(
        [int(c) for c in seq[-1].all_coeffs()]
    )


def main() -> None:
    A = associated_table(30)
    checked = 0
    for r in range(5, 31):
        F1, F2 = F_poly(A, r - 2), F_poly(A, r - 3)
        P, Q = U(F1), x * U(F2)
        w = W(P, Q)
        if w == 0:
            continue
        p = sp.Poly(w, x)
        nroots = real_roots(p)
        w0 = int(w.subs(x, 0))
        assert nroots == 0, (r, nroots)
        assert w0 < 0, (r, w0)
        checked += 1
    print(f"W(U F_{{r-2}}, x U F_{{r-3}}): OK ({checked} rows, r<=30)")


if __name__ == "__main__":
    main()
