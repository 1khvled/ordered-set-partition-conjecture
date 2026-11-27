#!/usr/bin/env python3
"""Guess short differential/lag recurrences and validate on held-out rows."""
import argparse
from itertools import product
import sympy as sp
from scan import ordered_rows


def deriv(c):
    return [k * c[k] for k in range(1, len(c))] or [0]


def mul_t(c, p=1):
    return [0] * p + list(c)


def add_scaled(dst, c, scale):
    if len(dst) < len(c):
        dst.extend([0] * (len(c) - len(dst)))
    for i, x in enumerate(c):
        dst[i] += scale * x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--j", type=int, default=3)
    ap.add_argument("--lag", type=int, default=3)
    ap.add_argument("--npow", type=int, default=2)
    ap.add_argument("--fit", type=int, default=24)
    ap.add_argument("--validate", type=int, default=45)
    args = ap.parse_args()
    rows = [c for _, c in ordered_rows(args.j, args.validate)]
    # operator = t^p D^q, with n^a multiplier; D = d/dt.
    basis = []
    for lag, a, p, q in product(range(args.lag + 1), range(args.npow + 1), range(3), range(3)):
        # Avoid identically zero high derivatives on the earliest rows later.
        basis.append((lag, a, p, q))
    equations = []
    for n in range(args.lag, args.fit + 1):
        nz = max(len(c) for c in rows[n-args.lag:n+1])
        for k in range(nz + 3):
            row = []
            for lag, a, p, q in basis:
                c = rows[n - lag]
                for _ in range(q):
                    c = deriv(c)
                c = mul_t(c, p)
                row.append((n ** a) * (c[k] if k < len(c) else 0))
            equations.append(row)
    M = sp.Matrix(equations)
    print("matrix", M.shape, "computing nullspace", flush=True)
    ns = M.nullspace()
    print("nullspace dimension", len(ns), flush=True)
    for idx, v in enumerate(ns[:10]):
        terms = [(basis[i], v[i]) for i in range(len(basis)) if v[i]]
        bad = []
        for n in range(args.fit + 1, args.validate + 1):
            acc = [0] * (len(rows[n]) + 10)
            for (lag, a, p, q), coef in terms:
                c = rows[n - lag]
                for _ in range(q):
                    c = deriv(c)
                c = mul_t(c, p)
                add_scaled(acc, c, int(coef) * n ** a)
            if any(acc):
                bad.append(n)
                if len(bad) >= 3:
                    break
        print("relation", idx, "terms", len(terms), "heldout_bad", bad)
        print(terms)


if __name__ == "__main__":
    main()
