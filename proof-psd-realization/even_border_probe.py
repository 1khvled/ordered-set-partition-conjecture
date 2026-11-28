"""Exact even-step PSD-bordering probe for j=2.

If M_{2m}=[[M_{2m-1},u],[u^T,alpha]], then
hat O_{2m}/hat O_{2m-1}=1+t alpha-t^2 u^T(I+tM_{2m-1})^{-1}u.
Thus alpha is the linear coefficient of the ratio and
s_m=(1+t alpha-ratio)/t^2 is a Stieltjes transform when the bordering exists.
"""
from __future__ import annotations

from fractions import Fraction

from explore import ordered_rows, det_tridiagonal
from kill_test import (
    lanczos_from_moments,
    moments_of_ratio,
    primitive,
    series_division,
    sub,
    trim,
)


def even_border_data(p_odd: list[int], p_even: list[int]):
    # R=hat O_even/hat O_odd = p_odd0*p_even/(p_even0*p_odd).
    num = [Fraction(p_odd[0] * x) for x in p_even]
    den = [Fraction(p_even[0] * x) for x in p_odd]
    coeffs = series_division(num, den, 3)
    alpha = coeffs[1]
    # S=(1+alpha*t-R)/t^2.
    # ((1+alpha*t)*den-num)/t^2, over den.
    product = [Fraction(0)] * (len(den) + 1)
    for i, x in enumerate(den):
        product[i] += x
        product[i + 1] += alpha * x
    raw = sub(product, num)
    if all(x == 0 for x in raw):
        return alpha, [Fraction(0)], den, [], [], Fraction(0)
    assert raw[0] == 0 and raw[1] == 0
    snum = trim(raw[2:])
    sden = den
    snum, sden = primitive(snum, sden)
    size = len(sden) - 1
    if size == 0:
        return alpha, snum, sden, [], [], Fraction(0)
    moments = moments_of_ratio(snum, sden, 2 * size + 2)
    a, b, last_norm = lanczos_from_moments(moments, size)
    assert last_norm == 0
    assert det_tridiagonal(a, b) == [Fraction(x, p_odd[0]) for x in p_odd]
    assert all(x > 0 for x in a) and all(x > 0 for x in b[1:size])
    return alpha, snum, sden, a, b, moments[0]


def main() -> None:
    rows = ordered_rows(2, 11)
    for m in range(1, 6):
        odd, even = rows[2 * m - 1], rows[2 * m]
        alpha, num, den, a, b, weight = even_border_data(odd, even)
        print(f"m={m} alpha={alpha}")
        print(f"  s_{m}(t)=({num})/({den})")
        if a:
            print(f"  T_diag={[str(x) for x in a]}")
            print(f"  T_beta2={[str(x) for x in b[1:len(a)]]}")
            print(f"  u_weight={weight}")
            beta_digits = max([len(str(x.denominator)) for x in b[1:len(a)]] or [0])
            print(f"  max denominator digits diag={max(len(str(x.denominator)) for x in a)} beta={beta_digits}")
        else:
            print("  empty bordering vector")


if __name__ == "__main__":
    main()
