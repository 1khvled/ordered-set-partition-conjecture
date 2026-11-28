"""Test whether degree-grouped O rows form a fixed/nested Jacobi chain.

For D_q(t)=det(I+tJ_q), a bordering gives
D_{q+1}=(1+a_q t)D_q-b_q t^2 D_{q-1}.
The probe solves for a_q,b_q from coefficients and checks all remaining ones.
"""
from __future__ import annotations

from fractions import Fraction

from explore import ordered_rows


def normalized(p: list[int]) -> list[Fraction]:
    return [Fraction(x, p[0]) for x in p]


def fit_step(p_prev: list[int], p: list[int], p_next: list[int]):
    """Fit D_next=(1+a t)D-b t^2 D_prev; return (a,b,residual)."""
    d0, d1, d2 = normalized(p_prev), normalized(p), normalized(p_next)
    # t^1: e1(next)=e1(cur)+a.
    a = d2[1] - d1[1]
    # t^2: e2(next)=e2(cur)+a*e1(cur)-b.
    b = (d1[2] if len(d1) > 2 else 0) + a * d1[1] - d2[2]
    out = [Fraction(0)] * len(d2)
    for k, x in enumerate(d1):
        out[k] += x
        if k + 1 < len(out):
            out[k + 1] += a * x
    for k, x in enumerate(d0):
        if k + 2 < len(out):
            out[k + 2] -= b * x
    residual = [x - y for x, y in zip(out, d2)]
    return a, b, residual


def main() -> None:
    failures = []
    passes = []
    for j in range(2, 13):
        rows = ordered_rows(j, 12 * j + j - 1)
        for s in range(j):
            chain = [rows[q * j + s] for q in range(12)]
            for q in range(1, 11):
                a, b, residual = fit_step(chain[q - 1], chain[q], chain[q + 1])
                if any(residual):
                    failures.append((j, s, q, a, b, residual))
                    break
            else:
                passes.append((j, s))
    print("passing residues:", passes)
    print("first failures:")
    for rec in failures[:30]:
        j, s, q, a, b, residual = rec
        print(f"  j={j} s={s} q={q} a={a} b={b} residual={residual}")
    print("residues tested:", sum(range(2, 13)), "passes:", len(passes), "failures:", len(failures))


if __name__ == "__main__":
    main()
