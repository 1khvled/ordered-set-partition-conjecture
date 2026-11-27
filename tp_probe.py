#!/usr/bin/env python3
"""Probe total positivity of the ordered-big-block coefficient matrix."""
import argparse
import json
import random
import time
from pathlib import Path
from flint import fmpz_mat
from scan import ordered_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jmin", type=int, default=3)
    ap.add_argument("--jmax", type=int, default=12)
    ap.add_argument("--nmax", type=int, default=120)
    ap.add_argument("--samples", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=108)
    ap.add_argument("--out", type=Path, default=Path("tp-probe.jsonl"))
    args = ap.parse_args()
    rng = random.Random(args.seed)
    start = time.monotonic()
    count = 0
    with args.out.open("a", encoding="utf-8") as out:
        for j in range(args.jmin, args.jmax + 1):
            rows = [c for _, c in ordered_rows(j, args.nmax)]
            def a(n, m):
                return rows[n][m] if m < len(rows[n]) else 0
            for _ in range(args.samples):
                k = rng.randint(2, min(8, args.nmax + 1, args.nmax // j + 1))
                rr = sorted(rng.sample(range(args.nmax + 1), k))
                cc = sorted(rng.sample(range(args.nmax // j + 1), k))
                det = int(fmpz_mat([[a(n, m) for m in cc] for n in rr]).det())
                count += 1
                if det < 0:
                    rec = {"j": j, "rows": rr, "cols": cc, "det": det}
                    out.write(json.dumps(rec) + "\n")
                    out.flush()
                    print(json.dumps(rec), flush=True)
            print(f"j={j}, minors={count}, elapsed={time.monotonic()-start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
