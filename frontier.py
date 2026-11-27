#!/usr/bin/env python3
"""Push past the public Issue #108 frontier with exact arithmetic.

New region:
  * j = 3,4,5: 141 <= n <= 600;
  * j = 6..60: 7j <= n <= 8j-1.

Every row gets an exact integer discriminant and certified Arb root
isolation.  A negative discriminant or a non-real root enclosure triggers
an exact Sturm recount.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from flint import ctx, fmpz_poly

from scan import newton_failure, nonreal_roots, ordered_rows, sturm_variations


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--small-nmin", type=int, default=141)
    ap.add_argument("--small-nmax", type=int, default=600)
    ap.add_argument("--jmin", type=int, default=6)
    ap.add_argument("--jmax", type=int, default=60)
    ap.add_argument("--prec", type=int, default=192)
    ap.add_argument("--out", type=Path, default=Path("frontier-results.jsonl"))
    args = ap.parse_args()
    ctx.prec = args.prec
    jobs = [(j, args.small_nmin, args.small_nmax) for j in (3, 4, 5)]
    jobs += [(j, 7 * j, 8 * j - 1) for j in range(args.jmin, args.jmax + 1)]
    done = set()
    if args.out.exists():
        for line in args.out.read_text(encoding="utf-8").splitlines():
            try:
                r = json.loads(line)
                done.add((r["j"], r["n"]))
            except Exception:
                pass
    start = time.monotonic()
    count = 0
    with args.out.open("a", encoding="utf-8") as out:
        for j, lo, hi in jobs:
            for n, coeffs in ordered_rows(j, hi):
                if n < lo or (j, n) in done:
                    continue
                p = fmpz_poly(coeffs)
                disc = p.discriminant() if p.degree() > 0 else 1
                rec = {
                    "j": j,
                    "n": n,
                    "degree": p.degree(),
                    "discriminant_sign": (disc > 0) - (disc < 0),
                    "discriminant_bits": abs(int(disc)).bit_length(),
                }
                nf = newton_failure(coeffs)
                if nf:
                    rec["newton_failure"] = nf
                bad = nonreal_roots(coeffs, args.prec)
                if bad or disc < 0 or nf:
                    vm = sturm_variations(coeffs, minus_inf=True)
                    v0 = sturm_variations(coeffs, 0)
                    vp = sturm_variations(coeffs)
                    rec["sturm"] = {
                        "minus_inf": vm,
                        "zero": v0,
                        "plus_inf": vp,
                        "negative_roots": vm - v0,
                        "nonpositive_distinct_roots": vm - vp,
                    }
                    rec["coefficients"] = coeffs
                if bad:
                    rec["nonreal"] = [
                        {"real": str(r.real), "imag": str(r.imag), "multiplicity": m}
                        for r, m in bad
                    ]
                rec["seconds"] = round(time.monotonic() - start, 3)
                out.write(json.dumps(rec) + "\n")
                out.flush()
                count += 1
                if bad or disc < 0 or nf:
                    print(json.dumps(rec), flush=True)
            print(f"finished j={j}, rows={count}, elapsed={time.monotonic()-start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
