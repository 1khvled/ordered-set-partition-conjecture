"""Exact obstruction certificate for proof-j2-rootchain.

No numerical root computation and no finite enumeration are used.
"""
from __future__ import annotations

import sympy as sp


t, x, a = sp.symbols("t x a", positive=True)


def U(f: sp.Expr) -> sp.Expr:
    return sp.expand(2 * f + t * sp.diff(f, t))


def Q4(av: sp.Expr) -> sp.Expr:
    return av + sp.Rational(13, 5) * av * (av + 1) * t


def Q5(av: sp.Expr) -> sp.Expr:
    return av + sp.Rational(103, 40) * av * (av + 1) * t


def Q6(av: sp.Expr) -> sp.Expr:
    return (
        av
        + sp.Rational(203, 20) * av * (av + 1) * t
        + 13 * av * (av + 1) * (av + 2) * t**2
    )


def main() -> None:
    print("V1: vertical recurrence (V), generic a")
    for name, q in (("Q4", Q4), ("Q5", Q5), ("Q6", Q6)):
        residual = sp.expand(a * q(a + 1) - (a + 1) * q(a) - t * sp.diff(q(a), t))
        print(f"  {name}: residual = {residual}")
        assert residual == 0

    print("V2: horizontal recurrence (M) at n=6, generic a")
    residual_m = sp.expand(
        Q6(a) - (sp.diff(t * Q5(a), t) + 5 * a * t * Q4(a + 1))
    )
    print(f"  residual = {residual_m}")
    assert residual_m == 0

    q = Q4(1)
    p = Q5(1)
    A = Q6(1)
    u = U(q)
    print("V3: base rows and coupled partner")
    print(f"  q = {q}")
    print(f"  p = {p}")
    print(f"  A = {A}")
    print(f"  Uq = {u}")
    assert sp.expand(A - (sp.diff(t * p, t) + 5 * t * u)) == 0

    R0 = sp.expand(A.subs(t, -x))
    C0 = sp.expand(sp.diff(x * sp.diff(R0, x), x))
    disc = sp.discriminant(R0, x)
    su = sp.Rational(5, 39)
    sp_root = sp.Rational(20, 103)
    zeta = sp.Rational(203, 3120)

    print("V4: exact root-location certificate for the coupled interlacings")
    print(f"  R0 = {R0}")
    print(f"  discriminant = {disc} > 0: {disc > 0}")
    print(f"  0 < su < sp: {0 < su < sp_root}")
    print(f"  R0(su) = {sp.factor(R0.subs(x, su))} < 0: {R0.subs(x, su) < 0}")
    print(f"  R0(sp) = {sp.factor(R0.subs(x, sp_root))} < 0: {R0.subs(x, sp_root) < 0}")
    assert disc > 0
    assert 0 < su < sp_root
    assert R0.subs(x, su) < 0
    assert R0.subs(x, sp_root) < 0

    print("V5: exact obstruction to rho_1 < zeta_1")
    print(f"  C0 = {C0}")
    print(f"  zeta_1 = {zeta}")
    print(f"  0 < zeta_1 < su: {0 < zeta < su}")
    rz = sp.factor(R0.subs(x, zeta))
    print(f"  R0(zeta_1) = {rz} > 0: {rz > 0}")
    assert sp.expand(C0.subs(x, zeta)) == 0
    assert 0 < zeta < su
    assert rz > 0

    # R0 has positive leading coefficient and two positive roots. The two
    # negative evaluations above put su,sp inside (rho_1,rho_2). Since
    # zeta<su<rho_2 and R0(zeta)>0, zeta must be left of rho_1.
    print("  conclusion: zeta_1 < rho_1 (exact sign certificate)")


if __name__ == "__main__":
    main()
