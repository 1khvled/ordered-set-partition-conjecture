#!/usr/bin/env python3
"""Exact scan for Alexandersson's ordered big-block conjecture.

Coefficients are Python integers.  Root decisions use FLINT's rigorous Arb
complex-root enclosures, never float64 roots.  A negative answer is checked
again by an exact Sturm sequence before it is reported.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

from flint import ctx, fmpq_poly, fmpz_poly


def ordered_rows(j: int, nmax: int):
    """Yield (n, coefficients) for O_{n,j}, coefficients in ascending order."""
    if j < 2:
        raise ValueError("j must be at least 2")
    rows: list[list[int]] = [[1]]
    yield 0, rows[0]
    for n in range(1, nmax + 1):
        deg = n // j
        row = [0] * (deg + 1)
        small = min(n, j - 1)
        # Choose the first block.  There are binom(n,r) ways to choose its
        # r labels; the remaining labels form an arbitrary ordered partition.
        # Blocks of size < j are unweighted, and larger blocks contribute t.
        for r in range(1, small + 1):
            c = math.comb(n, r)
            prev = rows[n - r]
            for m, a in enumerate(prev):
                row[m] += c * a
        for r in range(j, n + 1):
            c = math.comb(n, r)
            prev = rows[n - r]
            for m, a in enumerate(prev):
                row[m + 1] += c * a
        rows.append(row)
        yield n, row


def certified_complex_roots(coeffs: list[int], prec: int = 128):
    """Return Arb root balls, refining until each imaginary part is decided."""
    p = fmpz_poly(coeffs)
    if p.degree() <= 0:
        return []
    workprec = max(64, int(prec))
    for _ in range(8):
        with ctx.workprec(workprec):
            roots = p.complex_roots()
        undecided = [
            (r, mult) for r, mult in roots
            if not r.imag.is_zero()
            and not (r.imag.lower() > 0 or r.imag.upper() < 0)
        ]
        if not undecided:
            return roots
        workprec *= 2
    raise ArithmeticError("could not certify whether all roots are real")


def nonreal_roots(coeffs: list[int], prec: int = 128):
    """Return rigorously non-real root enclosures, or [] if all roots are real."""
    return [
        (r, mult) for r, mult in certified_complex_roots(coeffs, prec)
        if not r.imag.is_zero()
    ]


def sturm_variations(coeffs: list[int], x=None, *, minus_inf: bool = False):
    """Exact Sturm sign variations at x, +infinity (default), or -infinity."""
    p = fmpz_poly(coeffs)
    if p.degree() <= 0:
        return 0
    seq = [p, p.derivative()]
    while seq[-1].degree() > 0:
        # fmpz_poly's % is a pseudo-remainder over Z, so use Q[x] here.
        rem = -(fmpq_poly(seq[-2]) % fmpq_poly(seq[-1]))
        if rem == 0:
            break
        # Positive scaling preserves Sturm signs and controls coefficient size.
        rem = rem.numer()
        rem = rem // abs(rem.content())
        seq.append(rem)
    signs = []
    for q in seq:
        if minus_inf:
            v = int(q[q.degree()]) * (-1 if q.degree() % 2 else 1)
        elif x is None:
            v = int(q[q.degree()])
        else:
            v = int(q(x))
        if v:
            signs.append(1 if v > 0 else -1)
    return sum(a != b for a, b in zip(signs, signs[1:]))


def newton_failure(coeffs: list[int]):
    """Cheap necessary condition: ultra log-concavity with no internal zeros."""
    if not coeffs:
        return {"k": 0, "reason": "empty polynomial"}
    d = len(coeffs) - 1
    if d == 0:
        return None
    for k, a in enumerate(coeffs):
        if a <= 0:
            return {"k": k, "reason": "nonpositive coefficient"}
    for k in range(1, d):
        lhs = coeffs[k] ** 2 * k * (d - k)
        rhs = coeffs[k - 1] * coeffs[k + 1] * (k + 1) * (d - k + 1)
        if lhs < rhs:
            return {"k": k, "reason": "Newton ultra-log-concavity", "lhs": lhs, "rhs": rhs}
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jmin", type=int, default=2)
    ap.add_argument("--jmax", type=int, default=20)
    ap.add_argument("--nmax", type=int, default=150)
    ap.add_argument("--prec", type=int, default=128)
    ap.add_argument("--out", type=Path, default=Path("scan-results.jsonl"))
    ap.add_argument("--certify-all", action="store_true")
    args = ap.parse_args()
    ctx.prec = args.prec
    done = set()
    if args.out.exists():
        for line in args.out.read_text(encoding="utf-8").splitlines():
            try:
                rec = json.loads(line)
                done.add((rec["j"], rec["n"]))
            except Exception:
                pass
    start = time.monotonic()
    with args.out.open("a", encoding="utf-8") as out:
        for j in range(args.jmin, args.jmax + 1):
            for n, coeffs in ordered_rows(j, args.nmax):
                if (j, n) in done:
                    continue
                rec = {"j": j, "n": n, "degree": len(coeffs) - 1}
                nf = newton_failure(coeffs)
                if nf is not None:
                    rec["newton_failure"] = nf
                bad = nonreal_roots(coeffs, args.prec)
                if bad:
                    vminus = sturm_variations(coeffs, minus_inf=True)
                    vzero = sturm_variations(coeffs, 0)
                    vplus = sturm_variations(coeffs)
                    rec["nonreal"] = [
                        {"real": str(r.real), "imag": str(r.imag), "multiplicity": m}
                        for r, m in bad
                    ]
                    rec["sturm"] = {"at_minus_inf": vminus, "at_zero": vzero, "at_plus_inf": vplus}
                    rec["negative_roots"] = vminus - vzero
                    rec["nonpositive_roots"] = vminus - vplus
                    rec["coefficients"] = coeffs
                elif args.certify_all:
                    rec["sturm_nonpositive_roots"] = sturm_variations(coeffs, minus_inf=True) - sturm_variations(coeffs)
                rec["seconds"] = round(time.monotonic() - start, 3)
                out.write(json.dumps(rec) + "\n")
                out.flush()
                if bad or nf is not None:
                    print(json.dumps(rec), flush=True)
            print(f"finished j={j}, elapsed={time.monotonic()-start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
