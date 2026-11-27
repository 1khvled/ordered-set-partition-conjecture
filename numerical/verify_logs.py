#!/usr/bin/env python3
"""Audit frontier JSONL/GZip logs against an independent corrected recurrence."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

from flint import fmpz_poly

sys.set_int_max_str_digits(0)


def reference_rows(j: int, nmax: int):
    """The corrected two-sum recurrence, written independently of the scanner."""
    rows = [[1]]
    yield 0, rows[0]
    for n in range(1, nmax + 1):
        row = [0] * (n // j + 1)
        for r in range(1, min(n, j - 1) + 1):
            c = math.comb(n, r)
            for m, a in enumerate(rows[n - r]):
                row[m] += c * a
        for r in range(j, n + 1):
            c = math.comb(n, r)
            for m, a in enumerate(rows[n - r]):
                row[m + 1] += c * a
        rows.append(row)
        yield n, row


def read_jsonl(path: Path):
    if not path.exists():
        return
    with path.open(encoding="utf-8") as fp:
        for line_no, line in enumerate(fp, 1):
            if line.strip():
                try:
                    yield json.loads(line)
                except Exception as exc:
                    raise ValueError(f"{path}:{line_no}: {exc}") from exc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["small", "wide"], required=True)
    ap.add_argument("--jmin", type=int, required=True)
    ap.add_argument("--jmax", type=int, required=True)
    ap.add_argument("--nmin", type=int, default=141)
    ap.add_argument("--nmax", type=int, default=800)
    ap.add_argument("--wide-start-multiplier", type=int, default=7)
    ap.add_argument("--wide-end-multiplier", type=int, default=8)
    ap.add_argument("--rows", type=Path, required=True)
    ap.add_argument("--pairs", type=Path, required=True)
    ap.add_argument("--exact", type=Path, required=True)
    ap.add_argument("--failures", type=Path)
    ap.add_argument("--out", type=Path, default=Path("log-verification.json"))
    ap.add_argument("--allow-partial", action="store_true")
    ap.add_argument("--recompute-discriminants", type=int, default=0)
    args = ap.parse_args()
    errors = []
    row_records = list(read_jsonl(args.rows))
    pair_records = list(read_jsonl(args.pairs))
    failure_records = list(read_jsonl(args.failures)) if args.failures else []
    by_row = {}
    for rec in row_records:
        key = (rec["j"], rec["n"])
        if key in by_row:
            errors.append(f"duplicate row {key}")
        by_row[key] = rec
    by_pair = {}
    for rec in pair_records:
        key = (rec["j"], rec["left_n"], rec["right_n"])
        if key in by_pair:
            errors.append(f"duplicate pair {key}")
        by_pair[key] = rec
    expected_rows = set()
    expected_pairs = set()
    reference = {}
    for j in range(args.jmin, args.jmax + 1):
        lo, hi = ((args.nmin, args.nmax) if args.mode == "small"
                  else (args.wide_start_multiplier * j, args.wide_end_multiplier * j - 1))
        expected_rows.update((j, n) for n in range(lo, hi + 1))
        expected_pairs.update((j, n - 1, n) for n in range(lo, hi + 1))
        reference.update({(j, n): c for n, c in reference_rows(j, hi) if n >= lo - 1})
    missing_rows = sorted(expected_rows - set(by_row))
    extra_rows = sorted(set(by_row) - expected_rows)
    missing_pairs = sorted(expected_pairs - set(by_pair))
    extra_pairs = sorted(set(by_pair) - expected_pairs)
    if missing_rows and not args.allow_partial:
        errors.append(f"missing {len(missing_rows)} rows; first={missing_rows[:5]}")
    if missing_pairs and not args.allow_partial:
        errors.append(f"missing {len(missing_pairs)} pairs; first={missing_pairs[:5]}")
    if extra_rows:
        errors.append(f"unexpected rows: {extra_rows[:5]}")
    if extra_pairs:
        errors.append(f"unexpected pairs: {extra_pairs[:5]}")
    for key, rec in by_row.items():
        if not rec.get("ok"):
            errors.append(f"row not certified: {key}")
        if rec.get("discriminant_sign", 0) < 0:
            errors.append(f"negative discriminant: {key}")
        if not rec.get("newton", {}).get("ok"):
            errors.append(f"Newton failure: {key}")
        if rec.get("sturm_performed"):
            if not rec.get("sturm_ok"):
                errors.append(f"Sturm failure: {key}")
            st = rec.get("sturm", {})
            if not (st.get("real_distinct") == rec.get("degree")
                    and st.get("negative_distinct") == rec.get("degree")
                    and st.get("positive_distinct") == 0):
                errors.append(f"Sturm count details mismatch: {key}")
        arb = rec.get("arb", {})
        if not (arb.get("pairwise_disjoint") and arb.get("imag_contains_zero_all")
                and arb.get("all_real_by_disjoint_arb") and arb.get("root_order_certified")):
            errors.append(f"Arb all-real proof absent: {key}")
        if rec.get("degree") != len(reference.get(key, [])) - 1:
            errors.append(f"degree/recurrence mismatch: {key}")
    for key, rec in by_pair.items():
        if rec.get("ok") is not True:
            errors.append(f"pair not certified: {key}")
        d, e = rec["degree_f"], rec["degree_g"]
        expected_comparisons = 2 * d - 1 if e == d else 2 * d
        if rec.get("comparisons") != expected_comparisons:
            errors.append(f"comparison count mismatch: {key}")
        if rec.get("certified_comparisons") != expected_comparisons:
            errors.append(f"uncertified comparison: {key}")
        if rec.get("common_degree") != 0:
            errors.append(f"common roots require separate treatment: {key}")
    if failure_records:
        errors.append(f"{len(failure_records)} failure alerts present")
    exact_count = 0
    disc_recomputed = 0
    exact_keys = set()
    if args.exact.exists():
        with gzip.open(args.exact, "rt", encoding="utf-8") as fp:
            for line_no, line in enumerate(fp, 1):
                try:
                    rec = json.loads(line)
                except Exception as exc:
                    errors.append(f"exact log line {line_no}: {exc}")
                    continue
                key = (rec["j"], rec["n"])
                exact_keys.add(key)
                exact_count += 1
                if key not in reference:
                    errors.append(f"unexpected exact row {key}")
                    continue
                if rec["coefficients"] != reference[key]:
                    errors.append(f"exact coefficients differ from corrected recurrence: {key}")
                disc = int(rec["discriminant"])
                text = rec["discriminant"]
                if hashlib.sha256(text.encode()).hexdigest() != by_row.get(key, {}).get("discriminant_sha256"):
                    errors.append(f"discriminant hash mismatch: {key}")
                if abs(disc).bit_length() != by_row.get(key, {}).get("discriminant_bits"):
                    errors.append(f"discriminant bit count mismatch: {key}")
                if ((disc > 0) - (disc < 0)) != by_row.get(key, {}).get("discriminant_sign"):
                    errors.append(f"discriminant sign mismatch: {key}")
                if disc_recomputed < args.recompute_discriminants:
                    actual = int(fmpz_poly(rec["coefficients"]).discriminant())
                    if actual != disc:
                        errors.append(f"stored discriminant is wrong: {key}")
                    disc_recomputed += 1
    if exact_keys != set(by_row):
        errors.append(f"exact/row key mismatch: exact={len(exact_keys)}, rows={len(by_row)}")
    contiguous = {}
    for j in range(args.jmin, args.jmax + 1):
        lo, hi = ((args.nmin, args.nmax) if args.mode == "small"
                  else (args.wide_start_multiplier * j, args.wide_end_multiplier * j - 1))
        n = lo
        while (j, n) in by_row and (j, n - 1, n) in by_pair:
            n += 1
        contiguous[j] = {"first_n": lo, "last_contiguous_n": n - 1,
                         "complete": n == hi + 1}
    out = {"mode": args.mode, "jmin": args.jmin, "jmax": args.jmax,
           "wide_start_multiplier": args.wide_start_multiplier,
           "wide_end_multiplier": args.wide_end_multiplier,
           "rows": len(by_row), "pairs": len(by_pair), "exact_records": exact_count,
           "failure_alerts": len(failure_records), "discriminants_recomputed": disc_recomputed,
           "missing_rows": len(missing_rows), "missing_pairs": len(missing_pairs),
           "contiguous_ranges": contiguous, "errors": errors, "ok": not errors}
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out), flush=True)
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
