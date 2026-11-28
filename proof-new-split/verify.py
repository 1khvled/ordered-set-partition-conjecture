"""Exact checks for the A_{n,j} + t B_{n,j} split.

Finite evidence only; not a proof of the conjecture.
"""
from __future__ import annotations

import sympy as sp

t = sp.symbols("t")


def rows(j: int, nmax: int):
    import math
    R = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = sp.Integer(0)
        for r in range(1, min(n, j - 1) + 1):
            s += math.comb(n, r) * R[n - r]
        for r in range(j, n + 1):
            s += t * math.comb(n, r) * R[n - r]
        R.append(sp.expand(s))
    return R


def B_direct(R, n: int, j: int):
    if n < j:
        return None
    c = sp.binomial(n - 1, j - 1)
    return sp.expand(2 * sp.diff(R[n - 1], t) + c * (t - 1) * sp.diff(R[n - j], t))


def B_sum(R, n: int, j: int):
    import math
    s = sp.Integer(0)
    for r in range(j, n):
        s += math.comb(n - 1, r) * (R[n - r] + R[n - 1 - r])
    return sp.expand(s)


def W(f, g):
    return sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))


def nreal(poly: sp.Poly) -> int:
    seq = sp.sturm(poly)

    def sign(p, minus=False):
        sg = 1 if p.LC() > 0 else -1
        if minus and p.degree() % 2 == 1:
            sg = -sg
        return sg

    def var(minus):
        s = [sign(p, minus) for p in seq if not p.is_zero]
        return sum(a != b for a, b in zip(s, s[1:]))

    return var(True) - var(False)


def no_real_zeros(p) -> bool:
    q = sp.Poly(p, t)
    return (not q.is_zero) and nreal(q) == 0


def main() -> None:
    for j in range(2, 7):
        R = rows(j, 30)
        checked = 0
        for n in range(j, 31):
            lhs, rhs = B_direct(R, n, j), B_sum(R, n, j)
            assert sp.expand(lhs - rhs) == 0, (j, n)
            if n >= j + 1:
                assert all(c >= 0 for c in sp.Poly(rhs, t).all_coeffs()), (j, n, "B negative coeff")
            checked += 1
        print(f"j={j}: derivative-free identity and B>=0 OK ({checked} rows)")
    forall_a = 0
    for j in range(3, 7):
        R = rows(j, 25)
        for n in range(j + 1, 26):
            import math
            c = sp.binomial(n - 1, j - 1)
            A = n * R[n - 1]
            for r in range(2, j):
                A += math.comb(n, r) * R[n - r]
            A += c * t * R[n - j]
            w = W(R[n - 1], A)
            if no_real_zeros(w):
                forall_a += 1
            else:
                print("A interlace failure/sign change", j, n, int(w.subs(t, 0)))
    print(f"j>=3 O_{{n-1}} interlaces A_n: {forall_a} rows checked (see failures above)")
    R = rows(2, 25)
    ok = 0
    for n in range(6, 26):
        f = R[n - 1]
        g = sp.expand(t * B_direct(R, n, 2))
        w = W(f, g)
        assert no_real_zeros(w) and int(w.subs(t, 0)) < 0, (n, w.subs(t, 0))
        ok += 1
    print(f"j=2 O_{{n-1}} interlaces tB_n: OK ({ok} rows, n=6..25)")
    # obstruction at n=5 for the parameter-free differential split
    L5 = 5 * R[4] + 2 * t * sp.diff(R[4], t)
    DO3 = sp.diff((t - 1) * R[3], t)
    lam = sp.symbols("lam")
    p = sp.Poly(L5 + lam * t * DO3, t)
    disc = sp.factor(sp.discriminant(p, t))
    print("n=5 obstruction discriminant:", disc)


if __name__ == "__main__":
    main()
