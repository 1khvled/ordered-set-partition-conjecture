#!/usr/bin/env python3
"""Independent small-polynomial tests of the exact Sturm implementation."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import random

import sympy as sp
from flint import ctx, fmpz_poly

from certified_core import (SturmChain, check_real_rooted, interval_interlacing,
                            newton_check, ordered_rows, poly, root_balls,
                            squarefree_part, wronskian_crosscheck)

x = sp.Symbol("x")


def sympy_counts(coeffs: list[int]) -> dict[str, int]:
    p = sp.Poly(sum(int(a) * x ** i for i, a in enumerate(coeffs)), x)
    sq = sp.Poly(p.sqf_part(), x)
    real = int(sq.count_roots(-sp.oo, sp.oo))
    neg = int(sq.count_roots(-sp.oo, 0))
    pos = int(sq.count_roots(0, sp.oo))
    zero = int(p.eval(0) == 0)
    # SymPy's count_roots includes an endpoint root in both adjacent calls.
    return {"negative_distinct": neg - zero, "zero_roots": zero,
            "positive_distinct": pos - zero, "real_distinct": real,
            "squarefree_degree": sq.degree()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--random", type=int, default=250)
    ap.add_argument("--seed", type=int, default=108)
    ap.add_argument("--out", type=Path, default=Path("sturm-validation.json"))
    args = ap.parse_args()
    rng = random.Random(args.seed)
    cases = [
        [3], [1], [1, 1], [0, 0, 1], [1, 0, -2, 0, 1],
        [1, -2, 1], [-1, 0, 1], [1, 0, 1], [2, 3, 4],
        [1, -10, 35, -50, 24], [1, -6, 15, -20, 15, -6, 1],
        [1, -3, 3, -2, 1, -1],
    ]
    for _ in range(args.random):
        d = rng.randrange(0, 13)
        c = [rng.randrange(-8, 9) for _ in range(d + 1)]
        if c[-1] == 0:
            c[-1] = rng.choice([-7, -1, 1, 7])
        cases.append(c)
    mismatches = []
    method_mismatches = []
    negative_content_remainders = 0
    for coeffs in cases:
        p = poly(coeffs)
        chain = SturmChain.build(p, method="pseudo")
        rational_chain = SturmChain.build(p, method="rational")
        got = chain.summary()
        rational_got = rational_chain.summary()
        want = sympy_counts(coeffs)
        for key, val in want.items():
            if got[key] != val:
                mismatches.append({"coefficients": coeffs, "key": key,
                                   "got": got[key], "expected": val})
        if got != rational_got:
            method_mismatches.append({"coefficients": coeffs, "pseudo": got,
                                      "rational": rational_got})
        # Confirm every normalization used a positive divisor and preserved
        # the sign of the unnormalized rational remainder.
        if len(chain.polynomials) >= 3:
            from flint import fmpq_poly
            rem = -(fmpq_poly(chain.polynomials[-2]) % fmpq_poly(chain.polynomials[-1]))
            if rem.numer().content() < 0:
                negative_content_remainders += 1
    # End-to-end small-row and pair checks, including repeated/nonreal controls.
    precision_before = ctx.prec
    row_results = []
    for j in range(3, 9):
        prev_c = prev_cert = None
        for n, c in ordered_rows(j, 45):
            cert = check_real_rooted(c, prec=192, sturm=True)
            row_results.append({"j": j, "n": n, "ok": cert["ok"]})
            if not cert["ok"]:
                mismatches.append({"ordered_row": {"j": j, "n": n}, "certificate": cert})
            if prev_c is not None:
                pair = interval_interlacing(prev_c, c, prec=192, f_cert=prev_cert, g_cert=cert)
                if pair["ok"] is not True:
                    mismatches.append({"ordered_pair": {"j": j, "n": n - 1}, "result": pair})
            prev_c, prev_cert = c, cert
    wr = []
    wr_nonreal = []
    rows3 = dict(ordered_rows(3, 12))
    for n in range(0, 11, 2):
        rec = wronskian_crosscheck(rows3[n], rows3[n + 1], prec=192, sturm=True)
        wr.append({"n": n, "ok": rec["ok"]})
        if not rec["ok"]:
            # Deliberately not a primary-check failure: real-rooted Wronskians
            # are not being used as a universal interlacing criterion.
            wr_nonreal.append({"n": n, "result": rec})
    precision_after = ctx.prec
    known = [c for _, c in ordered_rows(3, 4)]
    expected = [[1], [1], [3], [12, 1], [66, 9]]
    rec = {"seed": args.seed, "random_cases": args.random, "cases": len(cases),
           "sturm_mismatches": len(mismatches), "mismatches": mismatches,
           "pseudo_rational_mismatches": len(method_mismatches),
           "method_mismatches": method_mismatches,
           "negative_content_remainders_exercised": negative_content_remainders,
           "known_j3_rows_ok": known == expected, "known_j3_rows": known,
           "ordered_rows_checked": len(row_results), "ordered_row_failures": sum(not r["ok"] for r in row_results),
           "global_precision_before": precision_before,
           "global_precision_after": precision_after,
           "global_precision_unchanged": precision_before == precision_after,
           "wronskian_crosschecks": wr,
           "wronskian_nonreal_crosschecks": wr_nonreal}
    args.out.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in rec.items() if k not in ("mismatches", "method_mismatches", "known_j3_rows", "wronskian_crosschecks", "wronskian_nonreal_crosschecks")}), flush=True)
    if mismatches or method_mismatches:
        print(json.dumps(mismatches[:3], indent=2), flush=True)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
