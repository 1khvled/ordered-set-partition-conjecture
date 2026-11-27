#!/usr/bin/env python3
"""Independent EGF and explicit ordered-partition checks of the recurrence.

For ordered set partitions with weight t on blocks of size at least j,

  F_j(x,t) = 1 / (1 - sum_{r=1}^{j-1} x^r/r!
                      - t (e^x - sum_{r=0}^{j-1} x^r/r!)).

Thus O_{n,j}(t) = n! [x^n] F_j(x,t).  All series arithmetic below is exact
Fraction arithmetic: it generically inverts the bivariate denominator and
does not call ordered_rows until the final compare.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path

from certified_core import ordered_rows


def egf_rows(j: int, nmax: int) -> list[list[int]]:
    fact = [math.factorial(i) for i in range(nmax + 1)]
    maxk = nmax // j
    # D = 1 - sum_{r<j} x^r/r! - t H, represented as D[n](t).
    den = [[Fraction(0) for _ in range(maxk + 1)] for _ in range(nmax + 1)]
    den[0][0] = 1
    for i in range(1, nmax + 1):
        if i < j:
            den[i][0] = Fraction(-1, fact[i])
        else:
            den[i][1] = Fraction(-1, fact[i])
    # Generic formal-series inverse D*F=1, truncated in x and t.
    inv = [[Fraction(0) for _ in range(maxk + 1)] for _ in range(nmax + 1)]
    inv[0][0] = 1
    for n in range(1, nmax + 1):
        for i in range(1, n + 1):
            for k, a in enumerate(den[i]):
                if not a:
                    continue
                for ell, b in enumerate(inv[n - i]):
                    if b and k + ell <= maxk:
                        inv[n][k + ell] -= a * b
    out = []
    for n in range(nmax + 1):
        coeffs = []
        for k in range(n // j + 1):
            v = fact[n] * inv[n][k]
            assert v.denominator == 1, (j, n, k, v)
            coeffs.append(v.numerator)
        out.append(coeffs)
    return out


def weight(partition: tuple[tuple[int, ...], ...], j: int) -> int:
    return sum(len(block) >= j for block in partition)


def explicit_ordered_partitions(j: int, nmax: int) -> list[Counter]:
    """Generate ordered set partitions by inserting label n in every way."""
    all_rows: list[Counter] = [Counter({0: 1})]
    partitions: list[tuple[tuple[int, ...], ...]] = [()]
    for n in range(1, nmax + 1):
        nxt = []
        counts: Counter = Counter()
        for partition in partitions:
            # Insert n into any existing block (append preserves sorted labels).
            for b, block in enumerate(partition):
                q = list(partition)
                q[b] = block + (n,)
                q = tuple(q)
                nxt.append(q)
                counts[weight(q, j)] += 1
            # Or create a singleton in any of the len(partition)+1 slots.
            for pos in range(len(partition) + 1):
                q = partition[:pos] + ((n,),) + partition[pos:]
                nxt.append(q)
                counts[weight(q, j)] += 1
        partitions = nxt
        all_rows.append(counts)
        if len(partitions) > 200_000:
            raise ValueError("explicit enumeration safety cap")
    return all_rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=40)
    ap.add_argument("--explicit-nmax", type=int, default=7)
    ap.add_argument("--jmin", type=int, default=3)
    ap.add_argument("--jmax", type=int, default=8)
    ap.add_argument("--out", type=Path, default=Path("egf-validation.jsonl"))
    args = ap.parse_args()
    with args.out.open("w", encoding="utf-8") as out:
        for j in range(args.jmin, args.jmax + 1):
            egf = egf_rows(j, args.nmax)
            rec_rows = [c for _, c in ordered_rows(j, args.nmax)]
            mismatches = [n for n in range(args.nmax + 1) if egf[n] != rec_rows[n]]
            exp = explicit_ordered_partitions(j, args.explicit_nmax)
            exp_rows = [[exp[n][k] for k in range(n // j + 1)] for n in range(args.explicit_nmax + 1)]
            exp_mismatch = [n for n in range(args.explicit_nmax + 1) if exp_rows[n] != egf[n]]
            rec = {"j": j, "egf_nmax": args.nmax, "explicit_nmax": args.explicit_nmax,
                   "egf_mismatches": mismatches, "explicit_mismatches": exp_mismatch,
                   "ok": not mismatches and not exp_mismatch,
                   "row_n0": egf[0], "row_n4": egf[4]}
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
