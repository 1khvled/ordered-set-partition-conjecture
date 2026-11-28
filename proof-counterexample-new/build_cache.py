#!/usr/bin/env python3
"""Generate and cache ordered-row coefficient lists for the counterexample hunt.

The recurrence in scan.ordered_rows builds rows from n=0, so reaching n=1200
for j=3 costs ~8 minutes.  Caching once makes hunt restarts cheap.
"""
from __future__ import annotations

import pickle
import sys
import time
from pathlib import Path

PARENT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PARENT))

from scan import ordered_rows

OUTDIR = Path(__file__).resolve().parent
CACHE = OUTDIR / "cache"
CACHE.mkdir(exist_ok=True)
LOG = OUTDIR / "cache_build.log"


def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(line + "\n")


def jobs():
    out = [(3, 1200), (4, 1200), (5, 1200)]
    for j in range(6, 61):
        out.append((j, 16 * j - 1))
    return out


def build(j, hi):
    path = CACHE / f"j{j}.pkl"
    if path.exists():
        try:
            with path.open("rb") as f:
                rows = pickle.load(f)
            if len(rows) > hi:
                log(f"j={j} cache hit rows={len(rows)} bytes={path.stat().st_size}")
                return
        except Exception as ex:
            log(f"j={j} cache invalid ({ex}); regenerating")
    t0 = time.monotonic()
    rows = []
    for n, coeffs in ordered_rows(j, hi):
        rows.append(coeffs)
        if n % 100 == 0:
            log(f"j={j} generating n={n}/{hi} elapsed={time.monotonic()-t0:.1f}s")
    tmp = path.with_suffix(".tmp")
    with tmp.open("wb") as f:
        pickle.dump(rows, f, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(path)
    log(f"j={j} cached rows={len(rows)} bytes={path.stat().st_size} elapsed={time.monotonic()-t0:.1f}s")


def main():
    for j, hi in jobs():
        build(j, hi)
    total = sum(p.stat().st_size for p in CACHE.glob("j*.pkl"))
    log(f"DONE cache_bytes={total} ({total/2**30:.3f} GiB)")


if __name__ == "__main__":
    main()
