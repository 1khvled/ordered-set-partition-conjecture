#!/usr/bin/env python3
"""Compact reproducible frontier audit; never stores full root enclosures."""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from flint import ctx, fmpz_poly
from scan import newton_failure, nonreal_roots, ordered_rows, sturm_variations
from interlace import check_pair


def digest(values):
    return hashlib.sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["j5", "j34", "wide"], required=True)
    ap.add_argument("--prec", type=int, default=192)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    args = ap.parse_args()
    if args.mode == "j5":
        jobs = [(5, 141, 800)]
    elif args.mode == "j34":
        jobs = [(3, 141, 800), (4, 141, 800)]
    else:
        jobs = [(j, 7*j, 12*j-1) for j in range(6, 61)]
    start = time.monotonic()
    rows = pairs = sturm_rows = failures = 0
    per_j = {}
    with args.out.open("w", encoding="utf-8") as out:
        for j, lo, hi in jobs:
            prev = None
            count = 0
            for n, coeffs in ordered_rows(j, hi):
                if n < lo:
                    prev = coeffs
                    continue
                p = fmpz_poly(coeffs)
                disc = p.discriminant() if p.degree() else 1
                rec = {
                    "kind": "row", "j": j, "n": n, "degree": p.degree(),
                    "sha256": digest(coeffs),
                    "discriminant_sign": (disc > 0) - (disc < 0),
                    "discriminant_sha256": hashlib.sha256(str(disc).encode()).hexdigest(),
                    "newton_ok": newton_failure(coeffs) is None,
                }
                bad = nonreal_roots(coeffs, args.prec)
                rec["arb_all_real"] = not bad
                if bad:
                    rec["nonreal"] = [{"real": str(r.real), "imag": str(r.imag), "multiplicity": m} for r, m in bad]
                # Exact Sturm on every wide row; on j=5 use a fixed schedule plus endpoints.
                if args.mode == "wide" or n == lo or n == hi or n % 15 == 0:
                    vm = sturm_variations(coeffs, minus_inf=True)
                    vp = sturm_variations(coeffs)
                    rec["sturm_distinct_real_roots"] = vm - vp
                    rec["sturm_complete"] = rec["sturm_distinct_real_roots"] == p.degree()
                    sturm_rows += 1
                if prev is not None:
                    ok, pair = check_pair(prev, coeffs, args.prec, wronskian=False)
                    pair.update({"kind": "pair", "j": j, "n": n - 1, "sha256_f": digest(prev), "sha256_g": digest(coeffs)})
                    out.write(json.dumps(pair, separators=(",", ":")) + "\n")
                    pairs += 1
                    if not ok:
                        failures += 1
                        pair["coefficients_f"] = prev
                        pair["coefficients_g"] = coeffs
                        out.write(json.dumps(pair, separators=(",", ":")) + "\n")
                row_bad = bad or disc < 0 or not rec["newton_ok"] or rec.get("sturm_complete") is False
                if row_bad:
                    failures += 1
                    rec["coefficients"] = coeffs
                out.write(json.dumps(rec, separators=(",", ":")) + "\n")
                rows += 1
                count += 1
                prev = coeffs
            per_j[str(j)] = count
            print(f"j={j} rows={count} total_rows={rows} failures={failures} elapsed={time.monotonic()-start:.1f}s", flush=True)
    summary = {
        "mode": args.mode, "precision_bits": args.prec, "rows": rows, "pairs": pairs,
        "exact_sturm_rows": sturm_rows, "failures": failures, "per_j": per_j,
        "seconds": round(time.monotonic()-start, 3), "ok": failures == 0,
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
