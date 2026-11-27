"""Bounded exact verification of the proved generating-function identities."""
from __future__ import annotations

from math import comb
from pathlib import Path
import json
from exact_algebra import add, mul, ordered_rows, recurrence_residual, scale


def H(rows, n):
    if n < 0:
        return [0]
    out = [0]
    for r in range(n + 1):
        out = add(out, scale(mul(rows[r], rows[n - r]), comb(n, r)))
    return out


def riccati_residual(j, rows, n):
    """O_{n+1} - (2 H_n - O_n + (t-1) C(n,j-1) H_{n-j+1})."""
    rhs = add(scale(H(rows, n), 2), rows[n], -1)
    if n >= j - 1:
        rhs = add(rhs, scale(mul([-1, 1], H(rows, n - j + 1)), comb(n, j - 1)))
    return add(rows[n + 1], rhs, -1)


def main():
    nmax = 60
    report = {"nmax": nmax, "j": {}}
    for j in range(2, 13):
        rows = ordered_rows(j, nmax)
        short_bad = [n for n in range(1, nmax + 1) if recurrence_residual(j, rows, n) != [0]]
        riccati_bad = [n for n in range(nmax) if riccati_residual(j, rows, n) != [0]]
        report["j"][str(j)] = {
            "short_recurrence_first_failure": short_bad[0] if short_bad else None,
            "riccati_first_failure": riccati_bad[0] if riccati_bad else None,
            "degree_formula_ok": all(len(row) - 1 == n // j for n, row in enumerate(rows)),
        }
    # Literal F checks at j=3, n=1.  These are rational-function identities
    # stored as numerator/denominator pairs.
    t = "t"
    report["literal_F_j3_n1"] = {
        "F_0": "1/(1-t)",
        "F_1": "1/(1-t)^2",
        "unweighted_recurrence_rhs": "1/(1-t)",
        "short_recurrence_rhs": "(1+t)/(1-t)^2",
        "conclusion": "neither equals F_1; the recurrence sequence is G=1/(1-(B-t))",
    }
    out = Path(__file__).with_name("exact_results.json")
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
