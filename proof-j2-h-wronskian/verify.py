"""Exact checks for proof-j2-h-wronskian.

Finite checks are finite evidence only, not an all-index proof.
"""
from __future__ import annotations
import sympy as sp

x, t, a = sp.symbols('x t a')

N = 30
A = [[0] * (N // 2 + 3) for _ in range(N + 1)]
A[0][0] = 1
for m in range(N):
    for r in range(len(A[0]) - 1):
        A[m + 1][r] += r * A[m][r]
        if r:
            A[m + 1][r] += m * A[m - 1][r - 1]


def rising(av, r):
    v = 1
    for i in range(r):
        v *= av + i
    return sp.Integer(v)


def Q(n, av, var=t):
    return sp.expand(sum(rising(av, r) * A[n][r] * var ** (r - 1)
                         for r in range(1, len(A[0])) if A[n][r]))


def Eop(f, var=t):
    return sp.expand(f + var * sp.diff(f, var))


def Top(f, var=t):
    return sp.expand(sp.diff(var * sp.diff(f, var), var))


def W(f, g, var=t):
    return sp.expand(sp.diff(f, var) * g - f * sp.diff(g, var))


def main():
    print('V1: identity H = E(L) - P*P''')
    P = sp.Function('P')(t)
    L = sp.expand(sp.diff(P, t)**2 - P*sp.diff(P, t, 2))
    H = W(P, Top(P), t)
    res = sp.expand(H - (Eop(L, t) - P*sp.diff(P, t, 2)))
    print('  residual =', res)
    assert res == 0

    print('V2: identity on actual Q_{n,a}, a=1,2, n<=12')
    for n in range(4, 13):
        for av in [1, 2]:
            Pn = Q(n, av)
            Lp = sp.expand(sp.diff(Pn, t)**2 - Pn*sp.diff(Pn, t, 2))
            H1 = W(Pn, Top(Pn), t)
            H2 = sp.expand(Eop(Lp, t) - Pn*sp.diff(Pn, t, 2))
            assert sp.expand(H1 - H2) == 0
    print('  residuals = 0')

    print('V3: direction sample n=6,a=1')
    P6 = Q(6, 1)
    H6 = W(P6, Top(P6), t)
    print('  H =', sp.factor(H6))
    assert sp.Poly(H6, t).count_roots(-sp.oo, sp.oo) == 0
    print('  no real roots, positive on R')

    print('V4: finite coefficient-positivity and real-root check (evidence only)')
    bad_c = []
    bad_r = []
    for n in range(4, 21):
        for av in [1, 2, 3]:
            Pn = Q(n, av)
            Hn = W(Pn, Top(Pn), t)
            coeffs = sp.Poly(Hn, t).all_coeffs()
            if any(c <= 0 for c in coeffs):
                bad_c.append((n, av))
            if sp.Poly(Hn, t).count_roots(-sp.oo, sp.oo):
                bad_r.append((n, av))
    print('  bad_coeff =', bad_c)
    print('  bad_real =', bad_r)
    assert not bad_c and not bad_r

    print('V5: boundary/degeneracy')
    for n in [4, 5]:
        Pn = Q(n, 1)
        d = sp.degree(Pn, t)
        print(f'  n={n}, deg P={d}, chain vacuous')
        assert d == 1

    print('V6: J-reduction obstruction from rootchain2 (reproduced)')
    J = (1 + t) * (1 + 2 * t) * (1 + 3 * t)
    Pj = Eop(J, t)
    Tj = Eop(sp.diff(Pj, t), t)
    Hj = W(Pj, Tj, t)
    val = Hj.subs(t, -sp.Rational(1, 10))
    print('  H(-1/10)=', val)
    assert val == -sp.Rational(2211, 625)

    print('ALL CHECKS PASSED (finite checks are not an all-index proof)')


if __name__ == '__main__':
    main()
