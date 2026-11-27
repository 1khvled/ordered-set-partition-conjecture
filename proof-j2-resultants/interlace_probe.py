"""Check: after removing the common factor t, the x-components
  Qt_nk(t) = (coeff of x^k in O_n(x,t)) / t
form an interlacing chain as k increases (for fixed n).
If so, any positive combination sum_k x^k Qt_nk is real-rooted with simple roots
by the interlacing-cone theorem, giving (D) for all x>0.
"""
from __future__ import annotations

import sympy as sp

x, t = sp.symbols("x t")
NMAX = 12


def rows_x(N):
    out = [sp.Integer(1)]
    for n in range(1, N + 1):
        p = n * x * out[n - 1]
        p += t * sum(sp.binomial(n, r) * out[n - r] for r in range(2, n + 1))
        out.append(sp.expand(p))
    return out


def real_roots(p):
    rr = sp.nroots(p, n=40, maxsteps=300)
    if any(abs(sp.im(r)) > sp.Float("1e-20") for r in rr):
        return None
    return sorted(float(sp.re(r)) for r in rr)


def interlace_strict(a, b, tol=1e-12):
    """a,b sorted lists of reals, strictly interlace (no common roots, alternate)."""
    merged = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    for i in range(len(merged) - 1):
        if abs(merged[i][0] - merged[i + 1][0]) < tol:
            return False  # shared root
        if merged[i][1] == merged[i + 1][1]:
            return False  # same source twice
    return True


def main():
    O = rows_x(NMAX)
    for n in range(4, NMAX + 1):
        Qs = []
        for k in range(0, n - 1):
            P = sp.expand(O[n].coeff(x, k))
            Q = sp.cancel(P / t)
            Q = sp.expand(Q)
            if sp.Poly(Q, t).degree() >= 1:
                Qs.append((k, sp.Poly(Q, t)))
        ok_all = True
        for (k1, p1), (k2, p2) in zip(Qs, Qs[1:]):
            r1, r2 = real_roots(p1.as_expr()), real_roots(p2.as_expr())
            if r1 is None or r2 is None:
                print(f"n={n} k={k1},{k2}: nonreal root", flush=True)
                ok_all = False
                continue
            ok = interlace_strict(r1, r2)
            if not ok:
                ok_all = False
                print(f"n={n} k={k1},{k2}: NOT interlacing. r1={r1} r2={r2}", flush=True)
        print(f"n={n}: chain interlacing = {ok_all} (len {len(Qs)})", flush=True)


if __name__ == "__main__":
    main()
