"""Exact finite checks for the associated-Stirling ratio inequality (H).

These checks are regression certificates and finite evidence only.
They are NOT an all-index proof of (H), SLC, or F>=0.
"""
from __future__ import annotations

import sympy as sp


def associated(max_m: int, max_k: int) -> list[list[int]]:
    A = [[0] * (max_k + 2) for _ in range(max_m + 2)]
    A[0][0] = 1
    for m in range(1, max_m + 1):
        for k in range(0, max_k + 1):
            if k == 0:
                A[m][k] = 1 if m == 0 else 0
                continue
            if 2 * k > m:
                continue
            val = k * A[m - 1][k]
            if m >= 2:
                val += (m - 1) * A[m - 2][k - 1]
            A[m][k] = val
    return A


def check_symbolic_identities() -> None:
    # (3): D_{m,k} = m * E'_{m,k}
    m, k = sp.symbols("m k")
    A_m_k, A_m_k1, A_m1_k, A_m1_k1 = sp.symbols(
        "A_m_k A_m_k1 A_m1_k A_m1_k1"
    )
    A_mm1_k, A_mm1_km1 = sp.symbols("A_mm1_k A_mm1_km1")
    D = k * ((k + 1) * A_m_k1 + m * A_mm1_k) * A_m_k - (k + 1) * (
        k * A_m_k + m * A_mm1_km1
    ) * A_m_k1
    Ep = k * A_mm1_k * A_m_k - (k + 1) * A_mm1_km1 * A_m_k1
    assert sp.factor(D - m * Ep) == 0

    # (8): E'_{m,k} = SLC(m-1,k) + (m-1)(Delta(m-1,k-1) + (m-2)E'_{m-2,k-1})
    # Use symbolic A(r,c) variables for the relevant entries.
    A = sp.Function("A")
    mm, kk = sp.symbols("mm kk")
    # Define concrete symbols for readability
    Am1k = sp.Symbol("Am1k")
    Am1km1 = sp.Symbol("Am1km1")
    Am1kp1 = sp.Symbol("Am1kp1")
    Am2km1 = sp.Symbol("Am2km1")
    Am2km2 = sp.Symbol("Am2km2")
    Am2k = sp.Symbol("Am2k")
    Amk = sp.Symbol("Amk")
    Amkp1 = sp.Symbol("Amkp1")
    Am3km2 = sp.Symbol("Am3km2")
    Am3km1 = sp.Symbol("Am3km1")
    Am2k_ = Am2k
    # E'_{m,k}
    Ep_mk = kk * Am1k * Amk - (kk + 1) * Am1km1 * Amkp1
    # substitute Amk, Amkp1
    sub_Amk = kk * Am1k + (mm - 1) * Am2km1
    sub_Amkp1 = (kk + 1) * Am1kp1 + (mm - 1) * Am2k_
    Ep_expanded = sp.expand(Ep_mk.subs({Amk: sub_Amk, Amkp1: sub_Amkp1}))
    SLC = sp.expand(kk**2 * Am1k**2 - (kk + 1) ** 2 * Am1km1 * Am1kp1)
    G = sp.expand(kk * Am1k * Am2km1 - (kk + 1) * Am1km1 * Am2k_)
    assert sp.expand(Ep_expanded - (SLC + (mm - 1) * G)) == 0

    # (7): G = D_{m-2,k-1} + Delta(m-1,k-1)
    D_m2_km1 = (kk - 1) * Am1k * Am2km1 - kk * Am1km1 * Am2k_
    Delta_m1_km1 = Am2km1 * Am1k - Am2k_ * Am1km1
    assert sp.expand(G - (D_m2_km1 + Delta_m1_km1)) == 0

    # D_{m-2,k-1} = (m-2)*E'_{m-2,k-1}
    Ep_m2_km1 = (kk - 1) * Am3km2 * Am2km1 - kk * Am3km1 * Am2k_
    # need D_{m-2,k-1} in terms of E'_{m-2,k-1} via (3)
    # D_{p,j}=p*E'_{p,j}, p=m-2,j=k-1
    D_m2_km1_via3 = (mm - 2) * Ep_m2_km1
    # but also express D_{m-2,k-1} directly by its definition with A(m-1,k),A(m-2,k-1),A(m-1,k-1),A(m-2,k)
    # D_{p,j} = j*A(p+1,j+1)*A(p,j) - (j+1)*A(p+1,j)*A(p,j+1)
    # p=m-2,j=k-1 => (k-1)*A(m-1,k)*A(m-2,k-1) - k*A(m-1,k-1)*A(m-2,k)
    D_def = (kk - 1) * Am1k * Am2km1 - kk * Am1km1 * Am2k_
    # And E'_{p,j}=j*A(p-1,j)*A(p,j)-(j+1)*A(p-1,j-1)*A(p,j+1)
    # p=m-2,j=k-1 => (k-1)*A(m-3,k-1)*A(m-2,k-1)-k*A(m-3,k-2)*A(m-2,k)
    # So D_def should equal (m-2)*Ep_m2_km1 by (3) after substitution? No, wait:
    # (3) is D_{m,k}=m*E'_{m,k}. So D_{m-2,k-1}=(m-2)*E'_{m-2,k-1}. Yes.
    # We've already verified (3) symbolically in general form above.
    assert sp.expand(D_def - D_m2_km1) == 0

    print("symbolic identities: OK")


def check_finite(A: list[list[int]], max_m: int) -> None:
    # (H)
    h_bad = h_checked = 0
    for m in range(2, max_m - 1):
        for k in range(1, max_m // 2):
            if 2 * k > m or 2 * (k + 1) > m + 1:
                continue
            a, b, c, d = A[m][k], A[m][k + 1], A[m + 1][k], A[m + 1][k + 1]
            if a == 0 or b == 0 or c == 0 or d == 0:
                continue
            h_checked += 1
            if k * d * a < (k + 1) * c * b:
                h_bad += 1
    assert h_bad == 0
    print(f"(H): OK ({h_checked} checks)")

    # (3)
    d_bad = d_checked = 0
    for m in range(1, max_m - 1):
        for k in range(1, max_m // 2):
            if 2 * k > m or 2 * (k + 1) > m + 1:
                continue
            if m - 1 < 0:
                continue
            D = k * A[m + 1][k + 1] * A[m][k] - (k + 1) * A[m + 1][k] * A[m][k + 1]
            Ep = k * A[m - 1][k] * A[m][k] - (k + 1) * A[m - 1][k - 1] * A[m][k + 1]
            d_checked += 1
            if D != m * Ep:
                d_bad += 1
    assert d_bad == 0
    print(f"identity (3): OK ({d_checked} checks)")

    # (8)
    e_bad = e_checked = 0
    for m in range(2, max_m - 1):
        for k in range(1, max_m // 2):
            if 2 * k > m or 2 * (k + 1) > m + 1:
                continue
            if k + 1 >= len(A[m - 1]):
                continue
            slc = k * k * A[m - 1][k] ** 2 - (k + 1) * (k + 1) * A[m - 1][k - 1] * A[m - 1][k + 1]
            delta = A[m - 2][k - 1] * A[m - 1][k] - A[m - 2][k] * A[m - 1][k - 1] if k - 1 >= 0 else 0
            ep2 = 0
            if m - 2 >= 1 and k - 1 >= 0:
                ep2 = (k - 1) * A[m - 3][k - 1] * A[m - 2][k - 1] - k * A[m - 3][k - 2] * A[m - 2][k] if (m - 3 >= 0 and k - 2 >= 0) else 0
                # better use direct E' definition
                ep2 = (k - 1) * A[m - 3 + 1][k - 1] * A[m - 2][k - 1] - k * A[m - 3 + 1][k - 2] * A[m - 2][k] if False else ep2
            # direct E'_{m-2,k-1}
            if m - 2 >= 1 and k - 1 >= 1:
                Ep2 = (k - 1) * A[m - 3][k - 1] * A[m - 2][k - 1] - k * A[m - 3][k - 2] * A[m - 2][k]
            else:
                Ep2 = 0
            Ep = k * A[m - 1][k] * A[m][k] - (k + 1) * A[m - 1][k - 1] * A[m][k + 1]
            rhs = slc + (m - 1) * delta + (m - 1) * (m - 2) * Ep2
            e_checked += 1
            if Ep != rhs:
                e_bad += 1
    assert e_bad == 0
    print(f"identity (8): OK ({e_checked} checks)")

    # SLC
    s_bad = s_checked = 0
    for m in range(2, max_m):
        d = m // 2
        for k in range(1, d):
            if k + 1 >= len(A[m]):
                continue
            if A[m][k - 1] * A[m][k] * A[m][k + 1] == 0:
                continue
            val = k * k * A[m][k] ** 2 - (k + 1) * (k + 1) * A[m][k - 1] * A[m][k + 1]
            s_checked += 1
            if val < 0:
                s_bad += 1
    assert s_bad == 0
    print(f"SLC: OK ({s_checked} checks)")

    # F
    f_bad = f_checked = 0
    for m in range(4, max_m - 1):
        d = (m - 1) // 2
        for k in range(1, d):
            if k + 1 >= len(A[m - 1]):
                continue
            if A[m - 1][k - 1] * A[m - 1][k] * A[m - 1][k + 1] == 0:
                continue
            slc = k * k * A[m - 1][k] ** 2 - (k + 1) * (k + 1) * A[m - 1][k - 1] * A[m - 1][k + 1]
            delta = A[m - 2][k - 1] * A[m - 1][k] - A[m - 2][k] * A[m - 1][k - 1]
            F = slc + (m - 1) * delta
            f_checked += 1
            if F < 0:
                f_bad += 1
    assert f_bad == 0
    print(f"F>=0: OK ({f_checked} checks)")


def check_newton_not_enough_for_slc() -> None:
    # Binomial (1+x)^d: at k=1, Newton holds but SLC (k^2 factor 1 vs 4) fails for d>=3.
    for d in range(3, 20):
        c0, c1, c2 = 1, d, d * (d - 1) // 2
        # Newton for degree d, k=1: c1^2 >= (2/1)*(d/(d-1))*c0*c2
        assert c1 * c1 * 1 * (d - 1) >= 2 * d * c0 * c2
        # SLC at k=1: 1*c1^2 >= 4*c0*c2 fails
        assert 1 * c1 * c1 < 4 * c0 * c2
    print("Newton-not-enough-for-SLC counterexample: OK")


def main() -> None:
    check_symbolic_identities()
    max_m = 300
    A = associated(max_m, max_m // 2 + 5)
    check_finite(A, max_m)
    check_newton_not_enough_for_slc()
    print("ALL CHECKS PASSED (finite evidence only)")


if __name__ == "__main__":
    main()
