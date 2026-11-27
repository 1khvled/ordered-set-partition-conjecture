#!/usr/bin/env python3
"""Produce precise per-j counts from validated row/pair JSONL logs."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path


def records(path: Path):
    if not path.exists():
        return
    with path.open(encoding="utf-8") as fp:
        for line in fp:
            if line.strip():
                yield json.loads(line)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=Path, required=True)
    ap.add_argument("--pairs", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("frontier-counts.json"))
    args = ap.parse_args()
    rows = defaultdict(list)
    pairs = defaultdict(list)
    for rec in records(args.rows):
        rows[rec["j"]].append(rec)
    for rec in records(args.pairs):
        pairs[rec["j"]].append(rec)
    out = {"row_file": str(args.rows), "pair_file": str(args.pairs), "j": {}}
    for j in sorted(set(rows) | set(pairs)):
        rr, pp = rows[j], pairs[j]
        d = {
            "rows": len(rr), "nmin": min((r["n"] for r in rr), default=None),
            "nmax": max((r["n"] for r in rr), default=None),
            "degree_min": min((r["degree"] for r in rr), default=None),
            "degree_max": max((r["degree"] for r in rr), default=None),
            "all_rows_certified": all(r["ok"] for r in rr),
            "all_discriminants_positive": all(r["discriminant_sign"] > 0 for r in rr),
            "discriminant_bits_min": min((r["discriminant_bits"] for r in rr), default=None),
            "discriminant_bits_max": max((r["discriminant_bits"] for r in rr), default=None),
            "newton_all_ok": all(r["newton"]["ok"] for r in rr),
            "newton_inequalities_checked": sum(r["newton"].get("checks", 0) for r in rr),
            "exact_sturm_rows": sum(bool(r["sturm_performed"]) for r in rr),
            "arb_only_rows": sum(not r["sturm_performed"] for r in rr),
            "pairs": len(pp),
            "left_nmin": min((r["left_n"] for r in pp), default=None),
            "right_nmax": max((r["right_n"] for r in pp), default=None),
            "all_pairs_certified": all(r["ok"] is True for r in pp),
            "pair_prec_bits_min": min((r["prec_bits"] for r in pp), default=None),
            "pair_prec_bits_max": max((r["prec_bits"] for r in pp), default=None),
            "interval_comparisons": sum(r["comparisons"] for r in pp),
            "wronskian_checks": sum("wronskian" in r for r in pp),
            "wronskian_nonreal": sum("wronskian" in r and not r["wronskian"]["ok"] for r in pp),
        }
        out["j"][str(j)] = d
    out["totals"] = {
        key: sum(d[key] for d in out["j"].values())
        for key in ("rows", "pairs", "newton_inequalities_checked", "exact_sturm_rows",
                    "arb_only_rows", "interval_comparisons", "wronskian_checks",
                    "wronskian_nonreal")
    }
    out["totals"]["all_rows_certified"] = all(d["all_rows_certified"] for d in out["j"].values())
    out["totals"]["all_pairs_certified"] = all(d["all_pairs_certified"] for d in out["j"].values())
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out["totals"]), flush=True)


if __name__ == "__main__":
    main()
