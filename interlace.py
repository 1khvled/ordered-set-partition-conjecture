#!/usr/bin/env python3
"""Certified consecutive-interlacing checks for O_{n,j}."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from flint import ctx, fmpz_poly

from scan import certified_complex_roots, ordered_rows


def root_enclosures(coeffs, prec):
    p = fmpz_poly(coeffs)
    roots = certified_complex_roots(coeffs, prec)
    if any(not r.imag.is_zero() for r, _ in roots):
        raise ValueError("polynomial is not real-rooted")
    vals = [r.real for r, mult in roots for _ in range(mult)]
    vals.sort(key=lambda a: float(a.mid()))
    # Sorting by midpoints is only a heuristic; certify the resulting order.
    if any(not certainly_leq(a, b) for a, b in zip(vals, vals[1:])):
        raise ArithmeticError("root order is not certified; increase precision")
    return p, vals


def certainly_leq(a, b):
    return bool(a.upper() <= b.lower())


def check_pair(coeffs_f, coeffs_g, prec=192, wronskian=False):
    f, a = root_enclosures(coeffs_f, prec)
    g, b = root_enclosures(coeffs_g, prec)
    d = len(a)
    e = len(b)
    if e not in (d, d + 1):
        return False, {"reason": "degree jump", "d": d, "e": e}
    if e == d:
        # f_1 <= g_1 <= f_2 <= ... <= f_d <= g_d.
        ok = all(certainly_leq(a[i], b[i]) for i in range(d))
        ok = ok and all(certainly_leq(b[i], a[i + 1]) for i in range(d - 1))
    else:
        # g_1 <= f_1 <= g_2 <= ... <= f_d <= g_{d+1}.
        ok = all(certainly_leq(b[i], a[i]) and certainly_leq(a[i], b[i + 1]) for i in range(d))
    rec = {"ok": bool(ok), "degree_f": d, "degree_g": e}
    if wronskian:
        w = f * g.derivative() - f.derivative() * g
        if w:
            roots = certified_complex_roots([int(c) for c in w], prec)
            bad = [(str(r.real), str(r.imag), m) for r, m in roots if not r.imag.is_zero()]
            rec["wronskian_degree"] = w.degree()
            rec["wronskian_all_real"] = not bad
            rec["wronskian_nonreal"] = bad
        else:
            rec["wronskian_all_real"] = True
        rec["common_root"] = f.gcd(g).degree() > 0
    return bool(ok), rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jmin", type=int, default=3)
    ap.add_argument("--jmax", type=int, default=20)
    ap.add_argument("--nmax", type=int, default=150)
    ap.add_argument("--prec", type=int, default=192)
    ap.add_argument("--wronskian", action="store_true")
    ap.add_argument("--out", type=Path, default=Path("interlace-results.jsonl"))
    args = ap.parse_args()
    ctx.prec = args.prec
    start = time.monotonic()
    failures = 0
    count = 0
    with args.out.open("a", encoding="utf-8") as out:
        for j in range(args.jmin, args.jmax + 1):
            prev = None
            for n, coeffs in ordered_rows(j, args.nmax):
                if prev is not None:
                    ok, rec = check_pair(prev, coeffs, args.prec, args.wronskian)
                    rec.update({"j": j, "n": n - 1, "seconds": round(time.monotonic() - start, 3)})
                    out.write(json.dumps(rec) + "\n")
                    out.flush()
                    count += 1
                    if not ok:
                        failures += 1
                        print(json.dumps(rec), flush=True)
                prev = coeffs
            print(f"finished j={j}, pairs={count}, failures={failures}, elapsed={time.monotonic()-start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
