"""Requested j=2 rank-one kill-test, using exact rational Lanczos only.

For r_m=(hat O_{2m+1}-hat O_{2m})/(t hat O_{2m}), a PSD rank-one update
would give r_m=v^T(I+tM)^{-1}v.  Its J-fraction is reconstructed from the
Taylor moments mu_k=(-1)^k[t^k]r_m, with no root computation.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm

from explore import ordered_rows, poly_inner, poly_sub, det_tridiagonal


def trim(p: list[Fraction]) -> list[Fraction]:
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def sub(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] -= x
    return trim(out)


def ratio_polynomials(p_even: list[int], p_odd: list[int]) -> tuple[list[Fraction], list[Fraction]]:
    """Return N,D with r=N/D, D(0) nonzero."""
    e0, o0 = p_even[0], p_odd[0]
    num = sub([Fraction(e0 * x) for x in p_odd], [Fraction(o0 * x) for x in p_even])
    assert num[0] == 0
    num = num[1:]  # divide by t
    den = [Fraction(o0 * x) for x in p_even]
    return trim(num), trim(den)


def primitive(num: list[Fraction], den: list[Fraction]) -> tuple[list[Fraction], list[Fraction]]:
    vals = num + den
    scale = 1
    for x in vals:
        scale = lcm(scale, x.denominator)
    ints = [int(x * scale) for x in vals]
    common = 0
    for x in ints:
        common = gcd(common, x)
    return [Fraction(x // common) for x in ints[:len(num)]], [Fraction(x // common) for x in ints[len(num):]]


def series_division(num: list[Fraction], den: list[Fraction], count: int) -> list[Fraction]:
    assert den[0]
    out: list[Fraction] = []
    for n in range(count):
        rhs = num[n] if n < len(num) else Fraction(0)
        rhs -= sum((den[i] * out[n - i] for i in range(1, min(n, len(den) - 1) + 1)), Fraction(0))
        out.append(rhs / den[0])
    return out


def moments_of_ratio(num: list[Fraction], den: list[Fraction], count: int) -> list[Fraction]:
    coeffs = series_division(num, den, count)
    return [(-1) ** k * c for k, c in enumerate(coeffs)]


def lanczos_from_moments(moments: list[Fraction], size: int) -> tuple[list[Fraction], list[Fraction], Fraction]:
    """Exact monic Lanczos; return alpha, beta (beta[0] unused), final norm."""
    v_prev: list[Fraction] = []
    v = [Fraction(1)]
    norm = poly_inner(v, v, moments)
    alpha: list[Fraction] = []
    beta = [Fraction(0)] * (size + 1)
    for n in range(size):
        xv = [Fraction(0)] + v
        a = poly_inner(xv, v, moments) / norm
        alpha.append(a)
        w = poly_sub(xv, v, a)
        if n:
            w = poly_sub(w, v_prev, beta[n])
        if n < size - 1:
            norm_next = poly_inner(w, w, moments)
            beta[n + 1] = norm_next / norm
            v_prev, v = v, w
            norm = norm_next
    return alpha, beta, poly_inner(w, w, moments)


def reconstruct_j_fraction(alpha: list[Fraction], beta: list[Fraction], weight: Fraction, order: int) -> list[Fraction]:
    """Taylor coefficients of weight/(1+a0 t-b1 t^2/(1+a1 t-...))."""
    # Build each tail denominator as a power series, from the inside out.
    tail = [Fraction(1)] + [Fraction(0)] * order
    for n in range(len(alpha) - 1, -1, -1):
        # tail currently represents 1/(1+a_{n+1}t-...); form
        # 1+a_n t-beta_{n+1}t^2*tail, then invert it.
        den = [Fraction(0)] * (order + 1)
        den[0] = Fraction(1)
        if order >= 1:
            den[1] = alpha[n]
        b = beta[n + 1] if n + 1 < len(beta) else Fraction(0)
        for k in range(order - 1):
            den[k + 2] -= b * tail[k]
        inv = series_division([Fraction(1)], den, order + 1)
        tail = inv
    return [weight * x for x in tail[:order]]


def main() -> None:
    rows = ordered_rows(2, 11)
    for m in range(1, 6):
        even, odd = rows[2 * m], rows[2 * m + 1]
        num, den = ratio_polynomials(even, odd)
        num, den = primitive(num, den)
        size = len(den) - 1
        moments = moments_of_ratio(num, den, 2 * size + 2)
        alpha, beta, last_norm = lanczos_from_moments(moments, size)
        assert last_norm == 0, (m, last_norm)
        assert det_tridiagonal(alpha, beta) == [Fraction(x, even[0]) for x in even]
        assert all(x > 0 for x in alpha) and all(x > 0 for x in beta[1:size])
        reconstructed = reconstruct_j_fraction(alpha, beta, moments[0], 2 * size + 2)
        expected_coeffs = [(-1) ** k * x for k, x in enumerate(moments_of_ratio(num, den, 2 * size + 2))]
        assert reconstructed == expected_coeffs
        print(f"m={m}")
        print(f"  r_{m}(t) = ({num})/({den})")
        print(f"  weight={moments[0]}")
        print(f"  alpha={[str(x) for x in alpha]}")
        print(f"  beta2={[str(x) for x in beta[1:size]]}")
        print(f"  max denominator digits alpha={max(len(str(x.denominator)) for x in alpha)} beta={max([len(str(x.denominator)) for x in beta[1:size]] or [0])}")


if __name__ == "__main__":
    main()
