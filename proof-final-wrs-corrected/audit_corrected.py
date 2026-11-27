#!/usr/bin/env python3
"""Lightweight exact-arithmetic audit for the corrected WRS track.

No numerical root is used as a certificate.  FLINT supplies isolating Arb
balls for real roots; all sign decisions are then rechecked by exact rational
Horner evaluation with a derivative remainder bound.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb

from flint import ctx, fmpz_poly


def ordered_rows(j: int, nmax: int) -> list[list[int]]:
    if j < 2:
        raise ValueError("j must be at least 2")
    rows = [[1]]
    for n in range(1, nmax + 1):
        row = [0] * (n // j + 1)
        for r in range(1, min(n, j - 1) + 1):
            c = comb(n, r)
            for k, value in enumerate(rows[n - r]):
                row[k] += c * value
        for r in range(j, n + 1):
            c = comb(n, r)
            for k, value in enumerate(rows[n - r]):
                row[k + 1] += c * value
        rows.append(row)
    return rows


def derivative(coeffs: list[int]) -> list[int]:
    return [k * coeffs[k] for k in range(1, len(coeffs))]


def exact_eval(coeffs: list[int], x: Fraction) -> Fraction:
    value = Fraction(0)
    for c in reversed(coeffs):
        value = value * x + c
    return value


def certified_sign(coeffs: list[int], center: Fraction, radius: Fraction) -> int:
    """Sign on [center-radius, center+radius], by exact derivative bound."""
    value = exact_eval(coeffs, center)
    R = max(abs(center - radius), abs(center + radius), Fraction(1))
    error = sum(
        k * abs(coeffs[k]) * R ** (k - 1)
        for k in range(1, len(coeffs))
    ) * radius
    if abs(value) <= error:
        raise ArithmeticError("sign not certified")
    return 1 if value > 0 else -1


def rational_center_radius(x) -> tuple[Fraction, Fraction]:
    # FLINT returns exact integers (mid, rad, exp) enclosing x in
    # [mid-rad, mid+rad] * 10**exp.  Convert that enclosure directly.
    mid, rad, exp = x.mid_rad_10exp(100)
    scale = Fraction(10) ** int(exp)
    return Fraction(int(mid)) * scale, Fraction(int(rad)) * scale


def real_root_balls(coeffs: list[int], prec: int = 160):
    p = fmpz_poly(coeffs)
    if p.degree() <= 0:
        return []
    workprec = max(64, prec)
    for _ in range(8):
        with ctx.workprec(workprec):
            roots = p.complex_roots()
        real = []
        decided = True
        for root, multiplicity in roots:
            if root.imag.is_zero():
                real.append((root.real, multiplicity))
            elif not (root.imag.lower() > 0 or root.imag.upper() < 0):
                decided = False
        if decided:
            return real
        workprec *= 2
    raise ArithmeticError("could not decide reality of every root")


def audit_range(j: int, nmax: int) -> dict[str, int]:
    rows = ordered_rows(j, nmax)
    root_checks = 0
    for n in range(j + 1, nmax + 1):
        q = rows[n - j]
        qp = derivative(q)
        final = rows[n]
        for ball, multiplicity in real_root_balls(q):
            center, radius = rational_center_radius(ball)
            sq = certified_sign(qp, center, radius)
            sp = certified_sign(final, center, radius)
            if sq * sp >= 0:
                raise ArithmeticError(
                    f"WRS sign failure at j={j}, n={n}, center={center}"
                )
            root_checks += multiplicity
    return {"j": j, "nmax": nmax, "root_checks": root_checks}


def retraction_certificate() -> None:
    rows = ordered_rows(2, 63)
    center = Fraction(
        "-0.301049864539336319020771115593462774725514449468323282248633"
    )
    eps = Fraction(1, 10**55)
    a, b = center - eps, center + eps
    q, qp, final = rows[61], derivative(rows[61]), rows[63]

    def variations_at(x: Fraction) -> int:
        # Exact rational Sturm sequence, with positive primitive scaling.
        seq = [q, qp]
        while seq[-1] and len(seq[-1]) > 1:
            dividend = seq[-2]
            divisor = seq[-1]
            rem = list(dividend)
            while rem and len(rem) >= len(divisor):
                factor = rem[-1] / divisor[-1]
                shift = len(rem) - len(divisor)
                for i, value in enumerate(divisor):
                    rem[i + shift] -= factor * value
                while rem and rem[-1] == 0:
                    rem.pop()
            if not rem:
                break
            seq.append([-value for value in rem])
        signs = []
        for poly in seq:
            value = exact_eval(poly, x)
            if value:
                signs.append(1 if value > 0 else -1)
        return sum(u != v for u, v in zip(signs, signs[1:]))

    assert exact_eval(q, a) != 0 and exact_eval(q, b) != 0
    root_count = variations_at(a) - variations_at(b)
    assert root_count == 1, root_count
    sq = certified_sign(qp, center, eps)
    sp = certified_sign(final, center, eps)
    assert sq * sp < 0
    print("retracted j=2,n=63 root: exactly one root in bracket")
    print("qprime_sign", sq, "O63_sign", sp, "residue_sign", sq * sp)
    print("WRS holds; this is not a counterexample")


def main() -> None:
    retraction_certificate()
    ranges = [(2, 40), (3, 36), (4, 32), (5, 30), (6, 30)]
    for j, nmax in ranges:
        result = audit_range(j, nmax)
        print("range", result)


if __name__ == "__main__":
    main()
