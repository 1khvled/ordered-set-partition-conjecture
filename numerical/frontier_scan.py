#!/usr/bin/env python3
"""Certified frontier scan, with exact discriminants and interval interlacing.

Primary real-rootedness certificate: rigorous FLINT/Arb complex-root balls,
pairwise disjoint, with every imaginary interval containing zero.  The
conjugate-overlap argument in certified_core proves all-real rootedness.
Exact Sturm recounts are also performed on a documented schedule.

Primary interlacing certificate: exact Sturm or the same disjoint-Arb proof
for both rows, followed by strict ordered Arb-ball inequalities. Wronskians
are sampled cross-checks only and never decide acceptance.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import time

# Exact discriminants at the frontier have far more than Python's default
# 4300 decimal digits. They are intentionally stored losslessly.
sys.set_int_max_str_digits(0)
from typing import Any

from flint import fmpz_poly

from certified_core import (SturmChain, check_real_rooted, interval_interlacing,
                            newton_check, ordered_rows, refine_nonreal_enclosures,
                            wronskian_crosscheck)


def json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items() if not k.startswith("_")}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)


def write_json_line(fp, rec: dict) -> None:
    fp.write(json.dumps(json_safe(rec), separators=(",", ":")) + "\n")
    fp.flush()


def load_done(path: Path, key: str) -> set:
    done = set()
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            try:
                rec = json.loads(line)
                done.add(rec[key])
            except Exception:
                pass
    return done


def exact_sturm_due(j: int, n: int, mode: str) -> bool:
    if mode == "wide":
        return True
    # All initial frontier rows, then a sparse deterministic high-degree net.
    return n <= 180 or (n % 100 == 0 and n <= 400)


def slim_cert(cert: dict) -> dict:
    return {k: json_safe(v) for k, v in cert.items() if not k.startswith("_")}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["small", "wide"], required=True)
    ap.add_argument("--jmin", type=int, default=3)
    ap.add_argument("--jmax", type=int, default=3)
    ap.add_argument("--nmin", type=int, default=141)
    ap.add_argument("--nmax", type=int, default=800)
    ap.add_argument("--wide-start-multiplier", type=int, default=7,
                    help="wide rows start at this multiple of j (default 7)")
    ap.add_argument("--wide-end-multiplier", type=int, default=8,
                    help="wide rows end at multiplier*j-1 (default 8j-1)")
    ap.add_argument("--prec", type=int, default=192)
    ap.add_argument("--out", type=Path, default=Path("frontier-results.jsonl"))
    ap.add_argument("--pairs-out", type=Path, default=Path("interlace-results.jsonl"))
    ap.add_argument("--exact-out", type=Path, default=Path("frontier-exact.jsonl.gz"))
    ap.add_argument("--failure-out", type=Path, default=Path("failure-alerts.jsonl"))
    ap.add_argument("--wronskian-every", type=int, default=20)
    args = ap.parse_args()
    if args.jmin < 3:
        raise SystemExit("j=2 is deliberately skipped")
    if args.mode == "wide" and args.wide_start_multiplier < 7:
        raise SystemExit("wide start must be at least 7j")
    if args.mode == "wide" and args.wide_end_multiplier <= args.wide_start_multiplier:
        raise SystemExit("wide endpoint multiplier must exceed the start multiplier")
    # All Arb operations use ctx.workprec internally; global precision is not
    # mutated by this scanner or by nested certificate calls.
    done_rows = load_done(args.out, "row_key")
    done_pairs = load_done(args.pairs_out, "pair_key")
    start = time.monotonic()
    stats = []
    with args.out.open("a", encoding="utf-8") as rows_out, \
         args.pairs_out.open("a", encoding="utf-8") as pairs_out, \
         gzip.open(args.exact_out, "at", encoding="utf-8") as exact_out, \
         args.failure_out.open("a", encoding="utf-8") as failure_out:
        for j in range(args.jmin, args.jmax + 1):
            if args.mode == "small":
                lo, hi = args.nmin, args.nmax
            else:
                lo, hi = args.wide_start_multiplier * j, args.wide_end_multiplier * j - 1
            job_start = time.monotonic()
            job = {"mode": args.mode, "j": j, "nmin": lo, "nmax": hi,
                   "rows_checked": 0, "rows_resumed": 0, "pairs_checked": 0,
                   "pairs_resumed": 0, "exact_sturm_rows": 0,
                   "arb_only_rows": 0, "wronskian_checks": 0,
                   "wronskian_nonreal": 0, "failures": 0,
                   "min_discriminant_sign": 1}
            prev_c = prev_cert = prev_data = prev_n = None
            # Generate once. The row at lo-1 is needed for the incoming pair.
            for n, coeffs in ordered_rows(j, hi):
                if n < lo - 1:
                    continue
                row_key = f"j{j}:n{n}"
                pair_key = f"j{j}:{n-1}->{n}"
                is_frontier = n >= lo
                due = exact_sturm_due(j, n, args.mode)
                method = "rational" if (args.mode == "wide" or len(coeffs) <= 22) else "pseudo"
                if not is_frontier:
                    # Support row for the incoming pair; not part of the claim.
                    cert = check_real_rooted(coeffs, args.prec, sturm=False,
                                             sturm_method=method, with_root_data=True)
                elif row_key in done_rows:
                    # Recompute quietly for pair continuity; do not duplicate logs.
                    cert = check_real_rooted(coeffs, args.prec, sturm=due,
                                             sturm_method=method, with_root_data=True)
                    job["rows_resumed"] += 1
                else:
                    cert = check_real_rooted(coeffs, args.prec, sturm=due,
                                             sturm_method=method, with_root_data=True)
                    newton = newton_check(coeffs)
                    disc = int(cert["_discriminant"])
                    disc_text = str(disc)
                    rec = {"j": j, "n": n, "row_key": row_key,
                           "degree": cert["degree"], "ok": cert["ok"],
                           "certification_method": cert["certification_method"],
                           "arb": cert["arb"], "arb_all_real_ok": cert["arb_all_real_ok"],
                           "sturm_performed": due,
                           "sturm_method": method if due else None,
                           "sturm": cert["sturm"], "sturm_ok": cert["sturm_ok"],
                           "discriminant_sign": (disc > 0) - (disc < 0),
                           "discriminant_bits": abs(disc).bit_length(),
                           "discriminant_sha256": hashlib.sha256(disc_text.encode()).hexdigest(),
                           "newton": newton, "seconds": round(time.monotonic() - start, 3)}
                    write_json_line(rows_out, rec)
                    write_json_line(exact_out, {"j": j, "n": n, "coefficients": coeffs,
                                                "discriminant": disc_text, "certificate": slim_cert(cert),
                                                "newton": newton})
                    job["rows_checked"] += 1
                    job["exact_sturm_rows"] += int(due)
                    job["arb_only_rows"] += int(not due)
                    job["min_discriminant_sign"] = min(job["min_discriminant_sign"], rec["discriminant_sign"])
                    if (not cert["ok"]) or (not newton["ok"]):
                        job["failures"] += 1
                        alert = {"kind": "row", "j": j, "n": n, "coefficients": coeffs,
                                 "discriminant": disc_text, "certificate": cert,
                                 "newton": newton}
                        alert["exact_sturm"] = SturmChain.build(fmpz_poly(coeffs), method="pseudo").summary()
                        alert["nonreal_enclosures"] = refine_nonreal_enclosures(coeffs, args.prec)
                        write_json_line(failure_out, alert)
                        print("FAILURE_ALERT " + json.dumps(json_safe(alert))[:2000], flush=True)
                data = cert.get("_root_data")
                if prev_c is not None and pair_key not in done_pairs:
                    pair = interval_interlacing(prev_c, coeffs, args.prec,
                                                f_cert=prev_cert, g_cert=cert,
                                                f_root_data=prev_data, g_root_data=data)
                    pair_rec = {"j": j, "left_n": n - 1, "right_n": n, "pair_key": pair_key,
                                "ok": pair["ok"], "reason": pair.get("reason"),
                                "degree_f": pair["degree_f"], "degree_g": pair["degree_g"],
                                "common_degree": pair["common_degree"],
                                "comparisons": pair.get("comparisons", 0),
                                "certified_comparisons": pair.get("certified", 0),
                                "prec_bits": pair.get("prec_bits", args.prec),
                                "f_real_rooted": pair.get("f_real_rooted"),
                                "g_real_rooted": pair.get("g_real_rooted"),
                                "seconds": round(time.monotonic() - start, 3)}
                    # Sampled Wronskians are recorded, including nonreal examples,
                    # but never alter pair acceptance.
                    left_index = (n - 1) - (lo - 1)
                    wd = pair["degree_f"] + pair["degree_g"] - 1
                    if ((left_index < 10 or left_index % args.wronskian_every == 0
                         or n == hi) and wd <= 80):
                        wr = wronskian_crosscheck(prev_c, coeffs, args.prec,
                                                  sturm=True)
                        pair_rec["wronskian"] = wr
                        job["wronskian_checks"] += 1
                        job["wronskian_nonreal"] += int(not wr["ok"])
                    write_json_line(pairs_out, pair_rec)
                    job["pairs_checked"] += 1
                    if pair["ok"] is not True:
                        job["failures"] += 1
                        alert = {"kind": "interlacing", "j": j, "left_n": n - 1,
                                 "right_n": n, "coefficients_left": prev_c,
                                 "coefficients_right": coeffs, "pair": pair,
                                 "left_exact_sturm": SturmChain.build(fmpz_poly(prev_c)).summary(),
                                 "right_exact_sturm": SturmChain.build(fmpz_poly(coeffs)).summary(),
                                 "left_nonreal": refine_nonreal_enclosures(prev_c, args.prec),
                                 "right_nonreal": refine_nonreal_enclosures(coeffs, args.prec)}
                        write_json_line(failure_out, alert)
                        print("FAILURE_ALERT " + json.dumps(json_safe(alert))[:2000], flush=True)
                elif prev_c is not None:
                    job["pairs_resumed"] += 1
                prev_c, prev_cert, prev_data, prev_n = coeffs, cert, data, n
                if is_frontier and (n < lo + 10 or n % 25 == 0 or n == hi):
                    print(json.dumps({"progress": True, **job, "n": n,
                                      "job_seconds": round(time.monotonic() - job_start, 1)}), flush=True)
            job["seconds"] = round(time.monotonic() - job_start, 3)
            stats.append(job)
            print("JOB_DONE " + json.dumps(job), flush=True)
    summary = {"mode": args.mode, "jmin": args.jmin, "jmax": args.jmax,
               "prec_bits": args.prec, "wide_start_multiplier": args.wide_start_multiplier,
               "wide_end_multiplier": args.wide_end_multiplier,
               "jobs": stats, "total_seconds": round(time.monotonic() - start, 3)}
    suffix = f"{args.mode}-j{args.jmin}" + (f"-{args.jmax}" if args.jmax != args.jmin else "")
    if args.mode == "wide" and (args.wide_start_multiplier, args.wide_end_multiplier) != (7, 8):
        suffix += f"-x{args.wide_start_multiplier}-{args.wide_end_multiplier}"
    Path(f"summary-{suffix}.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
