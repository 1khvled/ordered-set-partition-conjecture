"""Exact-root enclosure probe for the recurrence's compatibility window."""
from __future__ import annotations

from math import comb
import json
from flint import ctx, fmpz_poly
from exact_algebra import add, deriv, mul, ordered_rows, scale, t_deriv

ctx.prec = 160


def roots(p):
    p = trim_fmpz(p)
    if p.degree() == 0:
        return [0]
    rr = p.complex_roots()
    if any(not r.imag.is_zero() for r, _ in rr):
        return None
    vals = [r.real for r, mult in rr for _ in range(mult)]
    vals.sort(key=lambda x: float(x.mid()))
    return vals


def trim_fmpz(p):
    # fmpz_poly is already normalized.
    return p


def leq(a, b):
    return bool(a.upper() <= b.lower())


def lt(a, b):
    return bool(a.upper() < b.lower())


def pair_status(p, q):
    a, b = roots(p), roots(q)
    if a is None or b is None:
        return "nonreal"
    d, e = len(a) - 1, len(b) - 1
    if abs(d - e) > 1:
        return "degree_gap"
    # Check common roots exactly.
    if p.gcd(q).degree() > 0:
        return "common_root"
    if d == e:
        if d == 0:
            return "compatible"
        forward = all(lt(a[i], b[i]) and lt(b[i], a[i + 1]) for i in range(d - 1)) and lt(a[-1], b[-1])
        backward = all(lt(b[i], a[i]) and lt(a[i], b[i + 1]) for i in range(d - 1)) and lt(b[-1], a[-1])
        return "compatible" if forward or backward else "not_compatible"
    if d + 1 == e:
        return "compatible" if all(lt(b[i], a[i]) and lt(a[i], b[i + 1]) for i in range(d + 1)) else "not_compatible"
    return "compatible" if all(lt(a[i], b[i]) and lt(b[i], a[i + 1]) for i in range(e + 1)) else "not_compatible"


def main():
    failures = []
    summary = {}
    for j in range(2, 6):
        rows = ordered_rows(j, 18)
        counts = {"compatible": 0, "not_compatible": 0, "nonreal": 0, "degree_gap": 0, "common_root": 0}
        for n in range(j, 19):
            window = []
            for r in range(1, j + 1):
                if n - r >= 0:
                    window.append((f"O{n-r}", rows[n-r]))
            # R=n O_{n-1}+2T O_{n-1}
            R = add(scale(rows[n-1], n), scale(t_deriv(rows[n-1]), 2))
            window.append(("R", R))
            # Q=(t-1) D O_{n-j}=O+(t-1)O'
            Q = add(rows[n-j], mul([-1, 1], deriv(rows[n-j])))
            window.append(("Q", Q))
            window = [(name, p) for name, p in window if p != [0]]
            for i in range(len(window)):
                for k in range(i + 1, len(window)):
                    status = pair_status(fmpz_poly(window[i][1]), fmpz_poly(window[k][1]))
                    counts[status] += 1
                    if status != "compatible" and len(failures) < 50:
                        failures.append({"j": j, "n": n, "pair": [window[i][0], window[k][0]], "status": status})
        summary[str(j)] = counts
    report = {"summary": summary, "first_failures": failures}
    from pathlib import Path
    Path(__file__).with_name("window_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
