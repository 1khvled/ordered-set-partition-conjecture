"""Exact/rigorous checks for proof-associated-stirling-slc.

All finite checks are evidence only, not an all-index proof of SLC.
"""
from __future__ import annotations

import math
import mpmath
import sympy as sp


def associated(max_m: int, max_k: int) -> list[list[int]]:
    A = [[0] * (max_k + 2) for _ in range(max_m + 2)]
    A[0][0] = 1
    for m in range(1, max_m + 1):
        for k in range(0, max_k + 1):
            if k == 0:
                A[m][k] = 0
                continue
            if 2 * k > m:
                continue
            val = k * A[m - 1][k]
            if m >= 2:
                val += (m - 1) * A[m - 2][k - 1]
            A[m][k] = val
    return A


def check_slc(A: list[list[int]], max_m: int) -> None:
    bad = 0
    checked = 0
    worst = None
    for m in range(2, max_m + 1):
        d = m // 2
        for k in range(1, d):
            if A[m][k - 1] * A[m][k] * A[m][k + 1] == 0:
                continue
            slc = k * k * A[m][k] ** 2 - (k + 1) * (k + 1) * A[m][k - 1] * A[m][k + 1]
            checked += 1
            if slc < 0:
                bad += 1
            ratio = slc / ((k + 1) * (k + 1) * A[m][k - 1] * A[m][k + 1])
            if worst is None or ratio < worst[3]:
                worst = (m, k, slc, float(ratio))
    assert bad == 0
    print(f"SLC: OK ({checked} checks), worst ratio {worst}")


def check_gr_real_rooted(A: list[list[int]], max_m: int) -> None:
    mpmath.mp.dps = 80
    for r in range(2, max_m + 1):
        d = r // 2
        coeffs = [mpmath.mpf(math.factorial(k) ** 2 * A[r][k]) for k in range(d + 1)]
        roots = mpmath.polyroots(coeffs[::-1], maxsteps=300)
        nonreal = [rt for rt in roots if abs(mpmath.im(rt)) > mpmath.mpf("1e-30")]
        assert not nonreal, (r, len(nonreal))
    print(f"G_r real-rooted: OK through r={max_m}")


def check_c0_c2_control(A: list[list[int]], max_m: int) -> None:
    # SLC(m,k) = M^2*C2 + M*C1 + C0, M=m-1
    # C0 = k^2*SLC(m-1,k) + (k+1)^2*B*C  (weaker than SLC(m-1,k))
    # C2 = k^2*P^2 - (k+1)^2*Q*R is weaker than SLC(m-2,k-1)
    bad_c0 = bad_c2 = 0
    checked = 0
    for m in range(4, max_m + 1):
        M = m - 1
        for k in range(2, m // 2):
            if k + 1 > (m - 1) // 2:
                continue
            a = A[m - 1][k]
            b = A[m - 1][k - 1]
            c = A[m - 1][k + 1]
            p = A[m - 2][k - 1]
            q = A[m - 2][k - 2] if k - 2 >= 0 else 0
            r = A[m - 2][k]
            if a * b * c * p == 0:
                continue
            C0 = k**4 * a * a - (k * k - 1) * (k + 1) ** 2 * b * c
            slc_m1 = k * k * a * a - (k + 1) ** 2 * b * c
            assert C0 == k * k * slc_m1 + (k + 1) ** 2 * b * c
            C2 = k * k * p * p - (k + 1) ** 2 * q * r
            slc_m2_km1 = (k - 1) ** 2 * p * p - k * k * q * r
            checked += 1
            if slc_m1 >= 0 and C0 < 0:
                bad_c0 += 1
            if slc_m2_km1 >= 0 and C2 < 0:
                bad_c2 += 1
    assert bad_c0 == 0 and bad_c2 == 0
    print(f"C0/C2 control by SLC: OK ({checked} checks)")


def check_c1_positive(A: list[list[int]], max_m: int) -> None:
    bad = 0
    checked = 0
    worst = None
    for m in range(4, max_m + 1):
        for k in range(2, m // 2):
            if k + 1 > (m - 1) // 2:
                continue
            a = A[m - 1][k]
            b = A[m - 1][k - 1]
            c = A[m - 1][k + 1]
            p = A[m - 2][k - 1]
            q = A[m - 2][k - 2] if k - 2 >= 0 else 0
            r = A[m - 2][k]
            if a * b * c * p == 0:
                continue
            C1 = 2 * k**3 * a * p - (k + 1) ** 2 * ((k - 1) * b * r + (k + 1) * c * q)
            denom = (k + 1) ** 2 * ((k - 1) * b * r + (k + 1) * c * q)
            checked += 1
            if C1 < 0:
                bad += 1
            if denom > 0:
                ratio = C1 / denom
                if worst is None or ratio < worst[2]:
                    worst = (m, k, float(ratio))
    assert bad == 0
    print(f"C1>0: OK ({checked} checks), worst ratio {worst}")


def check_no_naive_lc_propagation() -> None:
    # Directly show the recurrence does not preserve LC by itself.
    b = [1, 5, 21, 21, 16, 12, 8, 5, 2]
    c = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    assert all(b[i] * b[i] >= b[i - 1] * b[i + 1] for i in range(1, len(b) - 1))
    assert all(c[i] * c[i] >= c[i - 1] * c[i + 1] for i in range(1, len(c) - 1))
    M = 1
    d = 8
    e = [0] * (d + 2)
    for k in range(1, d + 2):
        bk = b[k] if k < len(b) else 0
        ck1 = c[k - 1] if k - 1 < len(c) else 0
        e[k] = k * bk + M * k * k * ck1
    fails = [i for i in range(2, d) if e[i] * e[i] < e[i - 1] * e[i + 1]]
    assert fails
    print(f"naive LC propagation: counterexample confirmed at indices {fails}")


def check_slc_expansion_identity(A: list[list[int]], max_m: int) -> None:
    # Algebraic: SLC(m,k)=M^2*C2 + M*C1 + C0 exactly
    bad = 0
    checked = 0
    for m in range(4, max_m + 1):
        M = m - 1
        for k in range(2, m // 2):
            if k + 1 > (m - 1) // 2:
                continue
            a = A[m - 1][k]
            b = A[m - 1][k - 1]
            c_ = A[m - 1][k + 1]
            p = A[m - 2][k - 1]
            q = A[m - 2][k - 2] if k - 2 >= 0 else 0
            r = A[m - 2][k]
            if min(a, b, c_, p) == 0 and 2 * k > m:
                continue
            slc_m = k * k * A[m][k] ** 2 - (k + 1) ** 2 * A[m][k - 1] * A[m][k + 1]
            C2 = k * k * p * p - (k + 1) ** 2 * q * r
            C1 = 2 * k**3 * a * p - (k + 1) ** 2 * ((k - 1) * b * r + (k + 1) * c_ * q)
            C0 = k**4 * a * a - (k * k - 1) * (k + 1) ** 2 * b * c_
            checked += 1
            if slc_m != M * M * C2 + M * C1 + C0:
                bad += 1
    assert bad == 0
    print(f"SLC expansion identity: OK ({checked} checks)")


def main() -> None:
    max_m = 120
    A = associated(max_m, max_m // 2 + 5)
    check_slc(A, max_m)
    check_slc_expansion_identity(A, max_m)
    check_c0_c2_control(A, max_m)
    check_c1_positive(A, max_m)
    check_gr_real_rooted(A, 40)
    check_no_naive_lc_propagation()
    print("ALL CHECKS PASSED (finite evidence only)")


if __name__ == "__main__":
    main()
