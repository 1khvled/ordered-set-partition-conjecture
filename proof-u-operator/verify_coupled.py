"""Exact finite verification of the coupled G_r induction and (H)."""
from __future__ import annotations

import sympy as sp

x = sp.symbols("x")


def assoc(max_r: int):
    A = [[0] * (max_r + 2) for _ in range(max_r + 1)]
    A[0][0] = 1
    for r in range(2, max_r + 1):
        for k in range(1, r // 2 + 1):
            A[r][k] = k * A[r - 1][k] + (r - 1) * A[r - 2][k - 1]
    return A


def G(A, r: int):
    return sp.expand(sum(sp.factorial(k) * A[r][k] * x ** (k - 1)
                         for k in range(1, r // 2 + 1)))


def U(f):
    return sp.expand(2 * f + x * sp.diff(f, x))


def W(f, g):
    return sp.expand(sp.diff(f, x) * g - f * sp.diff(g, x))


def nreal(poly: sp.Poly) -> int:
    seq = sp.sturm(poly)
    def sgn(p, minus=False):
        sg = 1 if p.LC() > 0 else -1
        if minus and p.degree() % 2 == 1:
            sg = -sg
        return sg
    def var(minus):
        s = [sgn(p, minus) for p in seq if not p.is_zero]
        return sum(a != b for a, b in zip(s, s[1:]))
    return var(True) - var(False)


def no_real_zeros(poly) -> bool:
    p = sp.Poly(poly, x)
    return (not p.is_zero) and nreal(p) == 0


def main() -> None:
    A = assoc(31)
    Gs = {r: G(A, r) for r in range(3, 31)}
    for r in range(5, 31):
        f, g = Gs[r - 1], Gs[r]
        q = Gs[r - 2]
        # recurrence: G_r = (x G_{r-1})' + (r-1) x U G_{r-2}
        rec = sp.expand(sp.diff(x * f, x) + (r - 1) * x * U(q))
        assert sp.expand(g - rec) == 0, r
        assert nreal(sp.Poly(g, x)) == sp.Poly(g, x).degree(), (r, "real roots")
        w1 = W(f, g)
        assert no_real_zeros(w1) and int(w1.subs(x, 0)) < 0, (r, "G interlace", w1.subs(x, 0))
        w2 = W(U(f), g)
        assert no_real_zeros(w2) and int(w2.subs(x, 0)) < 0, (r, "UG interlace", w2.subs(x, 0))
    print("coupled G induction: OK for r=5..30")
    # (H): k A_{r,k+1} A_{r-1,k} >= (k+1) A_{r,k} A_{r-1,k+1}
    bad = 0
    for r in range(3, 31):
        for k in range(1, (r - 1) // 2 + 1):
            if k + 1 > (r // 2):
                continue
            lhs = k * A[r][k + 1] * A[r - 1][k]
            rhs = (k + 1) * A[r][k] * A[r - 1][k + 1]
            if lhs < rhs:
                bad += 1
    assert bad == 0
    print("(H): OK for r<=30")


if __name__ == "__main__":
    main()
