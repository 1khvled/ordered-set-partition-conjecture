#!/usr/bin/env python3
"""Independent exact rational root-interval validation of Arb interlacing."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from certified_core import ordered_rows

x = sp.Symbol("x")


def exact_intervals(coeffs: list[int], bits: int):
    p = sp.Poly(sum(int(a) * x ** i for i, a in enumerate(coeffs)), x)
    raw = p.intervals(eps=sp.Rational(1, 2 ** bits))
    vals = []
    for (a, b), mult in raw:
        vals.extend([(a, b)] * mult)
    vals.sort(key=lambda ab: ab[0])
    return p, vals


def leq_interval(a, b) -> bool:
    # Exact rational roots may be represented by degenerate intervals.
    return bool(a[1] <= b[0] or (a[0] == a[1] == b[0] == b[1]))


def validate_pair(cf: list[int], cg: list[int]) -> dict:
    d, e = len(cf) - 1, len(cg) - 1
    if e not in (d, d + 1):
        return {"ok": False, "reason": "degree jump"}
    for bits in (12, 24, 48, 96, 192):
        pf, a = exact_intervals(cf, bits)
        pg, b = exact_intervals(cg, bits)
        if len(a) != d or len(b) != e:
            return {"ok": False, "reason": "not all roots real", "bits": bits,
                    "roots_f": len(a), "roots_g": len(b)}
        if e == d:
            comparisons = ([(a[i], b[i]) for i in range(d)]
                           + [(b[i], a[i + 1]) for i in range(d - 1)])
        else:
            comparisons = ([(b[i], a[i]) for i in range(d)]
                           + [(a[i], b[i + 1]) for i in range(d)])
        if all(leq_interval(u, v) for u, v in comparisons):
            return {"ok": True, "bits": bits, "comparisons": len(comparisons),
                    "degree_f": d, "degree_g": e}
    return {"ok": None, "reason": "rational intervals did not separate", "bits": bits}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=0)
    ap.add_argument("--nmax", type=int, default=45)
    ap.add_argument("--jmin", type=int, default=3)
    ap.add_argument("--jmax", type=int, default=8)
    ap.add_argument("--out", type=Path, default=Path("rational-interlace-validation.json"))
    args = ap.parse_args()
    results = []
    failures = []
    for j in range(args.jmin, args.jmax + 1):
        prev = None
        count = 0
        for n, c in ordered_rows(j, args.nmax):
            if prev is not None and n - 1 >= args.nmin:
                rec = validate_pair(prev, c)
                rec.update({"j": j, "left_n": n - 1, "right_n": n})
                results.append(rec)
                count += 1
                if rec["ok"] is not True:
                    failures.append(rec)
            prev = c
        print(json.dumps({"j": j, "pairs": count, "failures": len(failures)}), flush=True)
    out = {"jmin": args.jmin, "jmax": args.jmax, "nmin": args.nmin, "nmax": args.nmax,
           "pairs": len(results), "failures": len(failures),
           "methods": "SymPy exact rational real-root isolation; independent of FLINT/Arb",
           "failure_records": failures}
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k != "failure_records"}), flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
