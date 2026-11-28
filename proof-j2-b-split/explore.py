"""Explore the j=2 B-split interlacing target.

Target: for n>=6, O_{n-1,2} interlaces t B_{n,2}.
Equivalent: at every negative root mu of O_{n-1},
  R_n(mu) := (n-1)*(1-mu)*O'_{n-2}(mu)/O'_{n-1}(mu) < 2.

Finite evidence only.
"""
from __future__ import annotations

import math
import sympy as sp

t = sp.symbols("t")


def rows(nmax: int):
    R = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = sp.Integer(0)
        for r in range(1, 2):  # j-1 = 1
            s += math.comb(n, r) * R[n - r]
        for r in range(2, n + 1):
            s += t * math.comb(n, r) * R[n - r]
        R.append(sp.expand(s))
    return R


def B_direct(R, n: int):
    # j=2, C(n-1,1)=n-1
    return sp.expand(2 * sp.diff(R[n - 1], t) + (n - 1) * (t - 1) * sp.diff(R[n - 2], t))


def W(f, g):
    return sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))


def neg_roots(poly_t):
    """Return sorted negative real roots as sympy RootOf/numeric."""
    p = sp.Poly(poly_t, t)
    rs = []
    for r in sp.real_roots(p):
        if r < 0:
            rs.append(r)
    return sorted(rs, key=lambda x: sp.N(x))


def main():
    N = 20
    R = rows(N)
    print("=== ratio R_n(mu) = (n-1)*(1-mu)*O'_{n-2}(mu)/O'_{n-1}(mu) at roots of O_{n-1} ===")
    for n in range(4, N + 1):
        On1 = R[n - 1]
        On2 = R[n - 2]
        dOn1 = sp.diff(On1, t)
        dOn2 = sp.diff(On2, t)
        roots = neg_roots(On1)
        if not roots:
            continue
        vals = []
        for mu in roots:
            val = sp.Rational(n - 1, 1) * (1 - mu) * dOn2.subs(t, mu) / dOn1.subs(t, mu)
            vals.append(sp.N(val, 20))
        print(f"n={n:2d} deg={int(sp.degree(On1,t)):2d} roots={len(roots)}  maxR={float(max(vals)):.12f}  vals={[float(v) for v in vals]}")

    print("\n=== check interlacing O_{n-1} vs t*B_n ===")
    for n in range(3, N + 1):
        f = R[n - 1]
        g = sp.expand(t * B_direct(R, n))
        w = W(f, g)
        wp = sp.Poly(w, t)
        # Sturm count
        seq = sp.sturm(wp)

        def sign_at(p, minus=False):
            lc = sp.LC(p)
            sg = 1 if lc > 0 else -1
            if minus and sp.degree(p) % 2 == 1:
                sg = -sg
            return sg

        def var(minus):
            s = [sign_at(p, minus) for p in seq if not p.is_zero]
            return sum(a != b for a, b in zip(s, s[1:]))

        nreal = var(True) - var(False)
        w0 = int(w.subs(t, 0))
        print(f"n={n:2d} nrealW={nreal} W(0)={w0}")

    print("\n=== ratio decomposition: O'_{n-1} recurrence? ===")
    # O_n = n O_{n-1} + t sum_{r=2}^{n} C(n,r) O_{n-r}
    # O'_n = n O'_{n-1} + sum_{r=2}^{n} C(n,r) O_{n-r} + t sum_{r=2}^{n} C(n,r) O'_{n-r}
    # Check identity: O'_n = n O'_{n-1} + (n*(n-1)/2)*? no, better compute directly
    for n in range(3, 12):
        On = R[n]
        On1 = R[n - 1]
        lhs = sp.diff(On, t)
        rhs = sp.Integer(0)
        for r in range(2, n + 1):
            rhs += math.comb(n, r) * (R[n - r] + t * sp.diff(R[n - r], t))
        rhs = sp.expand(n * sp.diff(On1, t) + rhs)
        print(n, sp.expand(lhs - rhs) == 0)


if __name__ == "__main__":
    main()
