#!/usr/bin/env python3
"""Killable compute worker for hunt.py (one JSON job per stdin line)."""
from __future__ import annotations

import argparse
import json
import sys
import traceback

import hunt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prec", type=int, default=64)
    ap.add_argument("--disc-degree-limit", type=int, default=64)
    args = ap.parse_args()
    hunt.ctx.prec = args.prec

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        job = {}
        try:
            job = json.loads(line)
            result = hunt.compute_job(job, args.prec, args.disc_degree_limit)
        except Exception as ex:
            traceback.print_exc()
            result = {"j": job.get("j"), "n": job.get("n"),
                      "error": f"{type(ex).__name__}: {ex}"}
        sys.stdout.write(json.dumps(result, separators=(",", ":"), default=str) + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
