"""Exact checks for the j=2 A+tB split interlacing target.

Finite evidence only; not an all-n proof.

Target (proof-new-split/NOTES.md §5):
  For all n>=6, O_{n-1,2} interlaces t B_{n,2}, where
    B_n = 2 O'_{n-1} + (n-1)(t-1) O'_{n-2}
        = sum_{r=2}^{n-1} C(n-1,r)(O_{n-r} + O_{n-1-r}).
  Equivalent: at every negative root mu of O_{n-1},
    R_n(mu) := (n-1)*(1-mu)*O'_{n-2}(mu)/O'_{n-1}(mu) < 2.

Run: py -3.14 verify.py
"""
from __future__ import annotations

import math
import sympy as sp

t = sp.symbols("t")


def rows(nmax: int):
    """j=2 rows O_n."""
    R = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = math.comb(n, 1) * R[n - 1]
        for r in range(2, n + 1):
            s += t * math.comb(n, r) * R[n - r]
        R.append(sp.expand(s))
    return R


def B_direct(R, n: int):
    return sp.expand(2 * sp.diff(R[n - 1], t) + (n - 1) * (t - 1) * sp.diff(R[n - 2], t))


def B_sum(R, n: int):
    s = sp.Integer(0)
    for r in range(2, n):
        s += math.comb(n - 1, r) * (R[n - r] + R[n - 1 - r])
    return sp.expand(s)


def W(f, g):
    return sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))


def no_real_zeros(p) -> bool:
    q = sp.Poly(sp.expand(p), t)
    if q.is_zero:
        return False
    seq = sp.sturm(q)

    def sign_at(poly, minus=False):
        lc = sp.LC(poly)
        sg = 1 if lc > 0 else -1
        if minus and sp.degree(poly) % 2 == 1:
            sg = -sg
        return sg

    def var(minus):
        s = [sign_at(poly, minus) for poly in seq if not poly.is_zero]
        return sum(a != b for a, b in zip(s, s[1:]))

    return (var(True) - var(False)) == 0


def neg_roots(poly_t):
    return sorted([r for r in sp.real_roots(sp.Poly(poly_t, t)) if r < 0], key=lambda x: sp.N(x))


def main() -> None:
    N = 25
    R = rows(N)

    # 1. B identity and coefficient positivity
    for n in range(2, N + 1):
        lhs, rhs = B_direct(R, n), B_sum(R, n)
        assert sp.expand(lhs - rhs) == 0, (n,)
        if n >= 3:
            assert all(c >= 0 for c in sp.Poly(rhs, t).all_coeffs()), (n,)
    print(f"1. B identity and B>=0 OK, n<= {N}")

    # 2. O_{n-1} interlaces tB_n, n=3..N
    for n in range(3, N + 1):
        f = R[n - 1]
        g = sp.expand(t * B_direct(R, n))
        w = W(f, g)
        assert no_real_zeros(w) and int(w.subs(t, 0)) < 0, (n,)
    print(f"2. O_{{n-1}} interlaces tB_n OK, n=3..{N}")

    # 3. equivalent B_n interlaces O_{n-1}
    for n in range(3, N + 1):
        w = W(B_direct(R, n), R[n - 1])
        assert no_real_zeros(w) and int(w.subs(t, 0)) < 0, (n,)
    print(f"3. B_n interlaces O_{{n-1}} OK, n=3..{N}")

    # 4. ratio R_n(mu) < 2 at roots of O_{n-1}
    for n in range(4, N + 1):
        On1, On2 = R[n - 1], R[n - 2]
        dOn1, dOn2 = sp.diff(On1, t), sp.diff(On2, t)
        for mu in neg_roots(On1):
            val = sp.N((n - 1) * (1 - mu) * dOn2.subs(t, mu) / dOn1.subs(t, mu))
            assert float(val) < 2.0, (n, mu, val)
    print(f"4. R_n(mu)<2 at roots of O_{{n-1}} OK, n=4..{N}")

    # 5. A_n interlacing follows from induction hypothesis (cone)
    # A_n = n O_{n-1} + (n-1) t O_{n-2}
    # If O_{n-2} <= O_{n-1}, then O_{n-1} <= t O_{n-2} in appropriate orientation
    # and O_{n-1} <= A_n. Check directly.
    for n in range(3, N + 1):
        A = sp.expand(n * R[n - 1] + (n - 1) * t * R[n - 2])
        w = W(R[n - 1], A)
        assert no_real_zeros(w), (n,)
    print(f"5. O_{{n-1}} interlaces A_n OK, n=3..{N}")

    # 6. individual B_n sum terms are NOT all positive at roots (cancellation needed)
    seen_neg = False
    for n in range(6, 12):
        On1 = R[n - 1]
        dOn1 = sp.diff(On1, t)
        for mu in neg_roots(On1):
            for r in range(2, n):
                term = math.comb(n - 1, r) * (R[n - r].subs(t, mu) + R[n - 1 - r].subs(t, mu)) / dOn1.subs(t, mu)
                if sp.N(term) < 0:
                    seen_neg = True
    assert seen_neg
    print("6. Individual B_n sum terms can be negative at roots of O_{n-1}: YES (cancellation required)")

    # 7. sharpness: R_n at largest root approaches 2
    vals = []
    for n in range(4, N + 1):
        On1, On2 = R[n - 1], R[n - 2]
        dOn1, dOn2 = sp.diff(On1, t), sp.diff(On2, t)
        roots = neg_roots(On1)
        mu = roots[-1]
        vals.append(float(sp.N((n - 1) * (1 - mu) * dOn2.subs(t, mu) / dOn1.subs(t, mu))))
    print(f"7. R_n at largest root, n=4..{N}: first={vals[0]:.6f} last={vals[-1]:.6f} increasing={all(a<b for a,b in zip(vals,vals[1:]))}")

    # 8. B_n interlaces B_{n-1}
    for n in range(4, N + 1):
        w = W(B_direct(R, n), B_direct(R, n - 1))
        assert no_real_zeros(w), (n,)
    print(f"8. B_n interlaces B_{{n-1}} OK, n=4..{N}")


if __name__ == "__main__":
    main()
