#!/usr/bin/env python3
"""Exact small-n checks for consecutive interlacing sign identities.

Checks:
- universal recurrence (R1)
- root evaluation identity (E1) / (Rho)
- first-block root identity (FBroot)
- NSC <=> interlacing (as in interlace.py convention)
- Wronskian root-sign equivalence (Wroot)
"""
from __future__ import annotations

import math
from fractions import Fraction

import sympy as sp

t = sp.Symbol('t')


def ordered_rows(j, nmax):
    rows = [[1]]
    yield 0, rows[0]
    for n in range(1, nmax + 1):
        deg = n // j
        row = [0] * (deg + 1)
        for r in range(1, min(n, j - 1) + 1):
            c = math.comb(n, r)
            prev = rows[n - r]
            for m, a in enumerate(prev):
                row[m] += c * a
        for r in range(j, n + 1):
            c = math.comb(n, r)
            prev = rows[n - r]
            for m, a in enumerate(prev):
                row[m + 1] += c * a
        rows.append(row)
        yield n, row


def poly_from_coeffs(coeffs):
    # coeffs are ascending: coeffs[i] * t^i
    return sp.Poly(sp.expand(sum(int(c) * t**i for i, c in enumerate(coeffs))), t)


def Tpoly(p):
    return sp.Poly(sp.expand(t * sp.diff(p.as_expr(), t)), t)


def Qop(p):
    return sp.Poly(sp.expand(p.as_expr() + (t - 1) * sp.diff(p.as_expr(), t)), t)


def Rn(n, p):
    return sp.Poly(sp.expand(n * p.as_expr() + 2 * Tpoly(p).as_expr()), t)


def Sn(n, j, O):
    expr = 0
    for r in range(2, j):
        if n - r >= 0:
            expr += math.comb(n, r) * O[n - r].as_expr()
    return sp.Poly(sp.expand(expr), t)


def real_roots(poly):
    # numeric roots, high precision; sufficient for sign/identity sanity checks
    if poly.degree() <= 0:
        return []
    roots = sp.nroots(poly.as_expr(), n=80)
    reals = []
    for r in roots:
        if abs(sp.im(r)) < sp.Rational(1,10**30):
            reals.append(sp.re(r))
    reals.sort(key=lambda x: float(sp.N(x,20)))
    return reals


def sgn(expr, prec=80):
    v = sp.N(expr, prec)
    # v may be complex with tiny imag
    if hasattr(v, 'as_real_imag'):
        rv, iv = v.as_real_imag()
        if abs(sp.N(iv,prec)) > sp.Rational(1,10**20):
            return None
        v = rv
    if abs(v) < sp.Rational(1,10**25):
        return 0
    return 1 if v > 0 else -1


def check():
    for j in range(2, 7):
        O = {}
        for n, coeffs in ordered_rows(j, 14):
            O[n] = poly_from_coeffs(coeffs)
        maxn = max(O)
        # (R1)
        for n in range(0, maxn):
            rhs = Rn(n + 1, O[n]).as_expr() + Sn(n + 1, j, O).as_expr()
            lag = n + 1 - j
            if lag >= 0:
                rhs += math.comb(n, j - 1) * sp.expand(t * Qop(O[lag]).as_expr())
            rhs = sp.Poly(sp.expand(rhs), t)
            assert rhs == O[n + 1], (j, n, 'R1')

        # root identities
        for n in range(0, maxn):
            P = O[n]
            Qp = O[n + 1]
            if P.degree() <= 0:
                continue
            roots = real_roots(P)
            assert len(roots) == P.degree(), (j, n, 'not simple/real?')
            Pexpr = P.as_expr(); Qexpr = Qp.as_expr()
            dP = sp.diff(Pexpr, t); dQ = sp.diff(Qexpr, t)
            for mu in roots:
                muN = sp.N(mu, 80)
                # E1
                e1_rhs = 2 * muN * dP.subs(t, muN) + Sn(n + 1, j, O).as_expr().subs(t, muN)
                lag = n + 1 - j
                if lag >= 0:
                    e1_rhs += math.comb(n, j - 1) * sp.expand(t * Qop(O[lag]).as_expr()).subs(t, muN)
                e1 = sp.N(Qexpr.subs(t, muN) - e1_rhs, 80)
                assert abs(e1) < 1e-30, (j, n, float(sp.N(mu,20)), e1, 'E1')

                # FBroot
                fb = sp.Integer(0)
                for r in range(2, j):
                    fb += math.comb(n + 1, r) * O[n + 1 - r].as_expr().subs(t, muN)
                for r in range(j, n + 2):
                    fb += muN * math.comb(n + 1, r) * O[n + 1 - r].as_expr().subs(t, muN)
                fb_err = sp.N(Qexpr.subs(t, muN) - fb, 80)
                assert abs(fb_err) < 1e-30, (j, n, 'FBroot')

                # Wroot
                Pp = dP.subs(t, muN); Qv = Qexpr.subs(t, muN)
                Wmu = -Pp * Qv
                lhs = Qv / Pp
                assert sgn(lhs) < 0, (j, n, mu, 'NSC sign')
                assert sgn(Wmu) > 0, (j, n, mu, 'Wroot')

        # NSC <=> interlacing via interlace.py ordering
        for n in range(0, maxn):
            P = O[n]
            Qp = O[n + 1]
            a = real_roots(P)
            b = real_roots(Qp)
            d, e = len(a), len(b)
            assert e in (d, d + 1)
            if e == d:
                inter = all(sp.simplify(a[i] - b[i]) <= 0 for i in range(d)) and all(sp.simplify(b[i] - a[i + 1]) <= 0 for i in range(d - 1))
            else:
                inter = all(sp.simplify(b[i] - a[i]) <= 0 and sp.simplify(a[i] - b[i + 1]) <= 0 for i in range(d))
            if d == 0:
                nsc = True
            else:
                nsc = all(sgn(Qp.as_expr().subs(t, mu) / sp.diff(P.as_expr(), t).subs(t, mu)) < 0 for mu in a)
            assert inter == nsc, (j, n, inter, nsc)
        print(f'j={j} ok')

if __name__ == '__main__':
    check()
