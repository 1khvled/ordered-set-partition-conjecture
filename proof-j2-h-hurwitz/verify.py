"""Exact checks for proof-j2-h-hurwitz.  Finite checks are evidence only."""
from __future__ import annotations
import sympy as sp
from common import Family, t

E, T, W = Family.E, Family.T, Family.W
F = Family(20)


def pos_on_R(p):
    """Exact certificate: no real roots, positive value at 0, positive leading coeff."""
    q = sp.Poly(sp.expand(p), t)
    return (q.count_roots(-sp.oo, sp.oo) == 0
            and q.eval(0) > 0 and sp.LC(q, t) > 0)


def main():
    print('V1 generic quadratic decomposition (symbolic residual)')
    c = sp.symbols('c')
    S = sp.Function('S')(t); K = sp.Function('K')(t)
    # work with expressions throughout
    Ee = lambda f: sp.expand(f + t * sp.diff(f, t))
    Uu = lambda f: sp.expand(2 * f + t * sp.diff(f, t))
    Ww = lambda f, g: sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))
    A = Ee(S); B = t * K
    alpha = Uu(Ee(sp.diff(S, t))); beta = Ee(Ee(K))
    H = sp.expand(Ww(A + c * B, alpha + c * beta))
    h0 = sp.expand(H.subs(c, 0))
    h2 = sp.expand(sp.diff(H, c, 2) / 2)
    h1 = sp.expand(sp.diff(H, c).subs(c, 0) - 2 * 0 * h2)
    res = sp.expand(H - (h0 + c * h1 + c**2 * h2))
    print('  residual =', res)
    assert res == 0
    print('  h0 =', h0)
    print('  h1 =', h1)
    print('  h2 =', h2)

    print('V2 identity h0 = W(E S, U T_S)')
    TSexpr = Ee(sp.diff(S, t))
    print('  residual =', sp.expand(h0 - Ww(A, Uu(TSexpr))))
    assert sp.expand(h0 - Ww(A, Uu(TSexpr))) == 0
    print('V3 identity h2 = W(tK, E^2 K)')
    print('  residual =', sp.expand(h2 - Ww(t * K, beta)))
    assert sp.expand(h2 - Ww(t * K, beta)) == 0

    print('V4 exact certificate for actual rows, 6<=n<=12, a in {1,2} (finite evidence)')
    bad = []
    for a in (1, 2):
        for n in range(6, 13):
            S = F.Qpoly(n - 1, a).as_expr()
            K = F.Qpoly(n - 2, a + 1).as_expr()
            Tprev = Ee(sp.diff(S, t))
            A = Ee(S); B = t * K
            alpha = Uu(Tprev); beta = Ee(Ee(K))
            h0 = Ww(A, alpha); h1 = Ww(A, beta) + Ww(B, alpha); h2 = Ww(B, beta)
            D = sp.expand(4 * h0 * h2 - h1**2)
            if not (pos_on_R(h0) and pos_on_R(h2) and pos_on_R(D)):
                bad.append((n, a))
    print('  failures =', bad)
    assert not bad

    print('V5 exact obstruction: coupled interlacings + h0,h2>0 are insufficient')
    q = sp.Rational
    S = (1 + t / q(8)) * (1 + t)
    T = 1 + t / q(7)
    K = (1 + t / q(9)) * (1 + t / q(3))
    print('  W(S,T) real roots =', sp.Poly(Ww(S, T), t).count_roots(-sp.oo, sp.oo))
    print('  W(K,S) real roots =', sp.Poly(Ww(K, S), t).count_roots(-sp.oo, sp.oo))
    assert sp.Poly(Ww(S, T), t).count_roots(-sp.oo, sp.oo) == 0
    assert sp.Poly(Ww(K, S), t).count_roots(-sp.oo, sp.oo) == 0
    A = Ee(S); B = t * K; alpha = Uu(T); beta = Ee(Ee(K))
    h0 = Ww(A, alpha); h2 = Ww(B, beta); h1 = Ww(A, beta) + Ww(B, alpha)
    D = sp.expand(4 * h0 * h2 - h1**2)
    print('  h0>0 on R?', pos_on_R(h0), 'h2>0 on R?', pos_on_R(h2))
    print('  D real roots =', sp.Poly(D, t).count_roots(-sp.oo, sp.oo))
    assert pos_on_R(h0) and pos_on_R(h2)
    assert sp.Poly(D, t).count_roots(-sp.oo, sp.oo) == 2
    Hc = sp.expand(h0 + h1 + h2)
    print('  c=1: H real roots =', sp.Poly(Hc, t).count_roots(-sp.oo, sp.oo))
    print('  H =', sp.factor(Hc))
    assert sp.Poly(Hc, t).count_roots(-sp.oo, sp.oo) == 2

    print('V6 boundary and direction samples')
    for n in (4, 5):
        P = F.Qpoly(n)
        print('  n=', n, 'deg P=', P.degree(), 'H=', F.H(n).as_expr())
    H6 = F.H(6)
    print('  H6=', sp.factor(H6.as_expr()), 'H6(0)=', H6.eval(0), '>0')
    assert H6.eval(0) > 0

    print('ALL EXACT CHECKS PASSED (finite evidence is not an all-index proof)')


if __name__ == '__main__':
    main()
