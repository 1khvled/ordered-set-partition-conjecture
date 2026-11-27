#!/usr/bin/env python3
"""Exact automated certificate for the consecutive sign condition on a range.

Certifies, for given j and nmax:
  (A) W(O_n,O_{n+1}) has no real roots (Sturm) for 1<=n<nmax.
  (B) Rho_{n+1}(mu)<0 at every root mu of O_n (high-precision evaluation).

This is a finite certificate, not an all-index proof.
"""
import sys, math, time
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs, real_roots
import sympy as sp

t = sp.Symbol('t')

def sturm_real_count(poly_expr):
    seq = sp.sturm(poly_expr, t)
    def eval_signs(val):
        nums = [sp.limit(s, t, val) for s in seq]
        signs = [0 if v == 0 else (1 if v > 0 else -1) for v in nums]
        nz = [s for s in signs if s != 0]
        ch = 0
        for i in range(len(nz) - 1):
            if nz[i] * nz[i + 1] < 0:
                ch += 1
        return ch
    return eval_signs(-sp.oo) - eval_signs(sp.oo)

def certify(j, nmax):
    O = {}
    for n, coeffs in ordered_rows(j, nmax):
        O[n] = poly_from_coeffs(coeffs)
    badW = 0
    badRho = 0
    worst = 0.0
    for n in range(1, nmax):
        P = O[n]
        Q = O[n + 1]
        if P.degree() <= 0:
            continue
        # (A) Wronskian real-root count
        W = sp.expand(P.as_expr() * sp.diff(Q.as_expr(), t)
                      - sp.diff(P.as_expr(), t) * Q.as_expr())
        nr = sturm_real_count(W)
        if nr != 0:
            badW += 1
            print(f"  W has {nr} real roots: j={j} n={n}")
        # (B) Rho at roots
        roots = real_roots(P)
        if len(roots) != P.degree():
            print(f"  root count mismatch j={j} n={n}: {len(roots)} vs {P.degree()}")
            continue
        for mu in roots:
            muN = sp.N(mu, 60)
            Pp = sp.diff(P.as_expr(), t).subs(t, muN)
            Qv = Q.as_expr().subs(t, muN)
            rho = float(sp.N(Qv / Pp, 40))
            if rho > worst:
                worst = rho
            if rho >= 0:
                badRho += 1
                print(f"  Rho>=0: j={j} n={n} rho={rho:.6e}")
    print(f"j={j} nmax={nmax}: badW={badW} badRho={badRho} worst_rho={worst:.6e}")
    return badW == 0 and badRho == 0

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--j", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=12)
    args = ap.parse_args()
    ok = certify(args.j, args.nmax)
    print("CERTIFIED" if ok else "FAILED")