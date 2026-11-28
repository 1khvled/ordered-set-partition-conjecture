"""Symbolic checks for the corrected C1 note.

This verifies valid identities and the exact obstruction to the failed optimization.
It does not assert an all-index C1 inequality and performs no finite enumeration.
"""
from __future__ import annotations

from fractions import Fraction

import sympy as sp


def main() -> None:
    n, k, y = sp.symbols("n k y", integer=True, positive=True)
    A, B, C, P, Q, R = sp.symbols("A B C P Q R", positive=True)

    slc_top = k**2 * A**2 - (k + 1) ** 2 * B * C
    slc_previous = (k - 1) ** 2 * P**2 - k**2 * Q * R

    # Valid cross-row ratio identity.
    rec_A = k * R + (n - 1) * P
    rec_B = (k - 1) * P + (n - 1) * Q
    cross_difference = sp.expand((k - 1) * rec_A * P - k * rec_B * R)
    expected_cross = sp.expand((n - 1) * ((k - 1) * P**2 - k * Q * R))
    assert sp.factor(cross_difference - expected_cross) == 0
    expected_plus_slack = sp.expand(
        (n - 1) * slc_previous / (k - 1)
        + (n - 1) * k * Q * R / (k - 1)
    )
    assert sp.factor(expected_cross - expected_plus_slack) == 0

    # Valid product constraint.
    product_gap = sp.expand(
        (k + 1) ** 2 * B * C * Q * R - (k - 1) ** 2 * A**2 * P**2
    )
    explicit_product_gap = sp.expand(-Q * R * slc_top - A**2 * slc_previous)
    assert sp.factor(product_gap - explicit_product_gap) == 0

    c = (k - 1) / (k + 1)
    r = (k - 1) / k
    alpha = (k + 1) ** 2 * (k - 1)
    beta = (k + 1) ** 3
    f = alpha * y + beta * c**2 / y

    # Correct monotonicity: f is decreasing, hence f(y) >= f(r) for y <= r.
    derivative_gap = sp.factor(beta * c**2 - alpha * r**2)
    assert sp.factor(derivative_gap - (k + 1) * (k - 1) ** 2 / k**2) == 0
    endpoint = sp.factor(alpha * r + beta * c**2 / r)
    assert sp.factor(endpoint - (2 * k**3 - 3 * k + 1 / k)) == 0
    difference = sp.factor(f - endpoint)
    expected_difference = sp.factor((r - y) * (beta * c**2 - alpha * r * y) / (y * r))
    assert sp.factor(difference - expected_difference) == 0

    # The constraints alone allow C1/X to tend to -infinity for every k>=2:
    # the singular coefficient is -beta*c^2/y with beta*c^2>0.
    obstruction = 2 * k**3 - alpha * y - beta * c**2 / y
    singular_coefficient = sp.factor(beta * c**2)
    assert sp.factor(singular_coefficient - (k + 1) * (k - 1) ** 2) == 0
    assert sp.limit(obstruction.subs(k, 2), y, 0, dir="+") == -sp.oo

    # Requested k=2 constraint example.
    y2 = Fraction(1, 10)
    z2 = Fraction(1, 30)
    assert y2 * z2 <= Fraction(1, 9)
    assert y2 <= Fraction(1, 2)
    c1_over_x = 16 - 9 * y2 - 27 * z2
    assert c1_over_x == Fraction(71, 5)
    assert c1_over_x > Fraction(11, 2)

    yneg = Fraction(1, 100)
    zneg = Fraction(100, 9)
    assert yneg * zneg == Fraction(1, 9)
    assert 16 - 9 * yneg - 27 * zneg == Fraction(-28409, 100)

    # Precise lower-root threshold for the product constraint alone.
    discriminant = sp.factor(k**6 - (k**2 - 1) ** 3)
    assert sp.factor(discriminant - (3 * k**4 - 3 * k**2 + 1)) == 0
    lower_root = (k**3 - sp.sqrt(discriminant)) / alpha
    root_polynomial = alpha * lower_root**2 - 2 * k**3 * lower_root + beta * c**2
    assert sp.simplify(root_polynomial) == 0

    # Valid Z=0 endpoint bound used in the partial cases.
    z0_bound = sp.factor(2 * k**3 - alpha * r)
    assert sp.factor(z0_bound - (k**3 + 2 * k - 1 / k)) == 0

    print("Corrected C1 symbolic checks: PASS")
    print("monotonicity corrected: f(y) >= f(r) on 0 < y <= r")
    print("formal constraint obstruction: C1/X -> -infinity as y -> 0+, z=c^2/y")
    print("k=2 example y=1/10, z=1/30 gives C1/X=71/5")
    print("negative feasible constraint example: y=1/100, z=100/9 gives C1/X=-28409/100")
    print("all-index C1 claim remains retracted")


if __name__ == "__main__":
    main()
