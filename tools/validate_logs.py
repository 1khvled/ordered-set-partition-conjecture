#!/usr/bin/env python3
"""Validate audit JSONL logs and emit a compact Markdown/JSON summary."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+", type=Path)
    ap.add_argument("--json-out", type=Path)
    ap.add_argument("--md-out", type=Path)
    args = ap.parse_args()
    by_j = defaultdict(lambda: {"count": 0, "min_n": None, "max_n": None, "max_degree": 0, "failures": []})
    total = 0
    for path in args.logs:
        if not path.exists():
            raise SystemExit(f"missing log: {path}")
        with path.open(encoding="utf-8") as f:
            for lineno, line in enumerate(f, 1):
                if not line.strip():
                    continue
                rec = json.loads(line)
                j, n = int(rec["j"]), int(rec["n"])
                d = by_j[j]
                d["count"] += 1
                d["min_n"] = n if d["min_n"] is None else min(d["min_n"], n)
                d["max_n"] = n if d["max_n"] is None else max(d["max_n"], n)
                d["max_degree"] = max(d["max_degree"], int(rec.get("degree", rec.get("degree_g", 0))))
                bad = []
                if rec.get("nonreal"):
                    bad.append("nonreal root enclosure")
                if rec.get("newton_failure"):
                    bad.append("Newton inequality")
                if rec.get("discriminant_sign", 1) < 0:
                    bad.append("negative discriminant")
                if rec.get("ok") is False:
                    bad.append("interlacing inequality")
                if rec.get("wronskian_all_real") is False:
                    bad.append("non-real Wronskian root (diagnostic)")
                if bad:
                    d["failures"].append({"file": str(path), "line": lineno, "n": n, "reasons": bad})
                total += 1
    summary = {
        "total_records": total,
        "by_j": {str(j): d for j, d in sorted(by_j.items())},
        "total_failure_records": sum(len(d["failures"]) for d in by_j.values()),
    }
    if args.json_out:
        args.json_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    lines = [
        "| j | records | n-range | max degree | flagged records |",
        "|---:|---:|:---|---:|---:|",
    ]
    for j, d in sorted(by_j.items()):
        lines.append(f"| {j} | {d['count']} | {d['min_n']}–{d['max_n']} | {d['max_degree']} | {len(d['failures'])} |")
    lines += ["", f"Total records: **{total}**", f"Flagged records: **{summary['total_failure_records']}**"]
    md = "\n".join(lines) + "\n"
    if args.md_out:
        args.md_out.write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
