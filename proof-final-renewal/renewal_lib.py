"""Exact-computation helpers for the threshold renewal study.

All ``O`` rows are computed with exact integer/symbolic arithmetic in ``t``
from the first-block convolution (no floating point).
"""
import sympy as sp
from math import comb

t = sp.symbols('t')


def EGF_zero_series(j, nmax):
    """Truncate e^z - sum_{r=0}^{j-1} z^r/r! to order nmax; return list c[n]."""
    c = [sp.Integer(0)] * (nmax + 1)
    for r in range(j, nmax + 1):
        c[r] = sp.Rational(1, sp.factorial(r))
    return c


def O_rows(j, nmax):
    """Exact rows O_{n,j}(t), n = 0..nmax, from first-block recurrence.

    O_0 = 1 ; O_n = sum_{r=1}^{j-1} binom(n,r) O_{n-r}
                    + t sum_{r=j}^{n} binom(n,r) O_{n-r}.
    """
    a = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        v = sp.Integer(0)
        for r in range(1, min(j - 1, n) + 1):
            v += comb(n, r) * a[n - r]
        for r in range(j, n + 1):
            v += comb(n, r) * t * a[n - r]
        a.append(sp.expand(v))
    return a


def binom_conv(A, B, N):
    """(A star B)_N = sum_{m=0}^N binom(N,m) A_m B_{N-m}."""
    out = sp.Integer(0)
    for m in range(N + 1):
        out += comb(N, m) * A[m] * B[N - m]
    return sp.expand(out)


def renewal_row(P, Q, n, j):
    """Q_n from (T4), using Q_m for m<n. Returns (RHS, correction)."""
    rhs = P[n]
    N = n - j + 1
    corr = sp.Integer(0)
    if N >= 0:
        corr = binom_conv(Q, P, N)
        rhs = sp.expand(rhs + (1 - t) * comb(n, j - 1) * corr)
    return rhs, corr


def verify_renewal(j, nmax):
    P = O_rows(j - 1, nmax)
    Q = O_rows(j, nmax)
    Q2 = [P[0]]
    for n in range(1, nmax + 1):
        rhs, _ = renewal_row(P, Q2, n, j)
        Q2.append(sp.expand(rhs))
    return all(sp.expand(Q2[n] - Q[n]) == 0 for n in range(nmax + 1))


def real_roots(poly):
    """All roots real? Return (ok, sorted float roots).

    Uses exact Sturm counting; numerical values only for ordering / display.
    """
    p = sp.Poly(sp.expand(poly), t)
    deg = p.degree()
    if deg <= 0:
        return True, []
    nreal = int(sp.polys.polytools.count_roots(p, -sp.oo, sp.oo))
    ok = (nreal == deg)
    roots = []
    try:
        roots = sorted(float(r) for r in sp.nroots(p, n=25, maxsteps=200))
    except Exception:
        # fall back to isolating intervals for real roots only
        try:
            ivs = sp.polys.polytools.intervals(p, eps=sp.Rational(1, 10**20))
            if isinstance(ivs, dict):
                roots = sorted(float(sp.Rational(a).evalf(30)) + 1e-12
                               for a, b in ivs.keys())
        except Exception:
            roots = []
    return ok, roots


def interlace(roots_a, roots_b):
    """a interlaces b (alternate, allowing the smaller-degree one)."""
    A, B = sorted(roots_a), sorted(roots_b)
    if len(A) > len(B):
        A, B = B, A
    if not A:
        return True
    # every root of A must lie in a distinct gap of B; for equal degree they alternate
    if len(A) == len(B):
        return all(A[i] <= B[i] <= A[i + 1] for i in range(len(A) - 1)) or \
               all(B[i] <= A[i] <= B[i + 1] for i in range(len(B) - 1))
    # |A| = |B|-1 case: A alternates in gaps of B
    if len(A) + 1 == len(B):
        return all(B[i] <= A[i] <= B[i + 1] for i in range(len(A)))
    return False


if __name__ == '__main__':
    for j in range(2, 7):
        print('j=%d renewal verified:' % j, verify_renewal(j, 15))
