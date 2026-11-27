"""Probe: O_n(x,t)=sum_k x^k P_n^{(k)}(t), with
P_n^{(k)}(t)=n![z^{n-k}](1-tC(z))^{-(k+1)}.
Check whether the polynomials P_n^{(k)} are real-rooted and interlace in k.
If so, any positive combination in x is real-rooted (interlacing cone).
"""
from __future__ import annotations

import sympy as sp

x, t = sp.symbols("x t")
NMAX = 10


def rows_x(N):
    out = [sp.Integer(1)]
    for n in range(1, N + 1):
        p = n * x * out[n - 1]
        p += t * sum(sp.binomial(n, r) * out[n - r] for r in range(2, n + 1))
        out.append(sp.expand(p))
    return out


def roots_of(p):
    roots = sp.nroots(p, n=30, maxsteps=200)
    real = all(abs(sp.im(r)) < sp.Float("1e-15") for r in roots)
    vals = sorted(sp.re(r) for r in roots)
    return real, vals


def main():
    O = rows_x(NMAX)
    for n in range(2, NMAX + 1):
        print(f"n={n}")
        prev = None
        for k in range(0, n + 1):
            P = sp.Poly(sp.expand(O[n].coeff(x, k)), t)
            d = P.degree()
            if d < 1:
                print(f"  k={k}: deg={d} coeffs={P.all_coeffs()}")
                continue
            real, vals = roots_of(P.as_expr())
            msg = ""
            if prev is not None and real and prev[0]:
                a, b = prev[1], vals
                # interlacing orientation check (roots of P vs roots of prev)
                msg = f" interlace={interlace(a, b)}"
            print(f"  k={k}: deg={d} real={real} roots={[round(float(v),4) for v in vals]}{msg}")
            prev = (real, vals)


def interlace(a, b):
    # a,b sorted lists; check b interlaces a (alternating)
    m, n = len(a), len(b)
    if abs(m - n) > 1:
        return False
    merged = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    # no ties allowed
    for i in range(len(merged) - 1):
        if merged[i][0] == merged[i + 1][0]:
            return False
    # strict alternation of source labels
    return all(merged[i][1] != merged[i + 1][1] for i in range(len(merged) - 1))


if __name__ == "__main__":
    main()
