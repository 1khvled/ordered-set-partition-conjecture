"""Exact residuals for the all-index proof of (Sign).

The first D2C premise is NOT proved here.  No numerical root checks are used
as proof.  Finite Sturm checks are printed only as labelled finite evidence
for the coupled theorem and for the still-open first premise.
"""
import sympy as sp


t, x = sp.symbols("t x", real=True)


def astirling_table(mmax):
    rmax = mmax // 2 + 2
    A = [[0] * rmax for _ in range(mmax + 1)]
    A[0][0] = 1
    for m in range(mmax):
        for r in range(rmax - 1):
            A[m + 1][r] += r * A[m][r]
            if r:
                A[m + 1][r] += m * A[m - 1][r - 1]
    return A


def rising(a, r):
    out = 1
    for k in range(r):
        out *= a + k
    return out


def Q(A, m, a, z):
    out = 0
    for r in range(1, len(A[0])):
        if A[m][r]:
            out += rising(a, r) * A[m][r] * z ** (r - 1)
    return sp.expand(out)


def B(f, g):
    """Opposite-oriented Wronskian in x: f*g' - f'*g."""
    return sp.expand(f * sp.diff(g, x) - sp.diff(f, x) * g)


def sep(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main():
    A = astirling_table(20)

    sep("V1: associated-Stirling recurrence (M) and vertical identity (V)")
    bad_m = bad_v = 0
    checked = 0
    for m in range(4, 16):
        for a in range(1, 6):
            q = Q(A, m, a, t)
            qm1 = Q(A, m - 1, a, t)
            qm2 = Q(A, m - 2, a + 1, t)
            qa1 = Q(A, m, a + 1, t)
            res_m = sp.expand(q - (sp.diff(t * qm1, t) + (m - 1) * a * t * qm2))
            res_v = sp.expand(qa1 - q - sp.diff(t * q, t) / a)
            bad_m += int(res_m != 0)
            bad_v += int(res_v != 0)
            checked += 1
    print(f"checked={checked}  M_failures={bad_m}  V_failures={bad_v}")

    sep("V2: generic sum-of-squares residual (SOS)")
    a, n = sp.symbols("a n", positive=True)
    R0, R1, R2 = sp.symbols("R0 R1 R2")
    V0 = ((a + 1) * R0 + x * R1) / a
    V1 = ((a + 2) * R1 + x * R2) / a
    Bvr = V0 * R1 - V1 * R0
    G = sp.expand(x * Bvr + (n + 1) * x * V0**2 + V0**2 - R0 * V0)
    M = 1 + (n + 1) * x * (a + 1)
    rhs = (
        a * x**2 * (R1**2 - R0 * R2)
        + (a + 1) * M * (R0 + x * R1 / (a + 1)) ** 2
        + a * x**2 * R1**2 / (a + 1)
    )
    sos_res = sp.factor(a**2 * G - rhs)
    print("generic SOS residual =", sos_res)

    sep("V3: generic refined residue identity (3.4)")
    H0, H1 = sp.symbols("H0 H1")
    S0 = sp.expand(a * (V0 - R0 - n * x * H0))
    S1 = sp.expand(a * (V1 - R1 - n * (H0 + x * H1)))
    Bsv = S0 * V1 - S1 * V0
    Bvh = V0 * H1 - V1 * H0
    E = sp.expand(x * Bsv + S0 * V0 + (n + 1) * a * x * V0**2)
    refined_rhs = sp.expand(a * (G + n * x**2 * Bvh))
    print("generic refined residual =", sp.factor(E - refined_rhs))

    sep("V4: root reduction P*V = -E modulo U")
    P0 = sp.expand(S0 + x * S1 - (n + 1) * a * x * V0)
    # At a root of U=2V+xV', use xV'=-2V.  Since V1 is displayed in
    # terms of R2, reduce modulo U as a polynomial in R2.
    U_generic = sp.expand(2 * V0 + x * V1)
    root_res = sp.factor(sp.rem(sp.expand(P0 * V0 + E), U_generic, R2))
    print("generic root residual =", root_res)

    sep("V5: exact actual-polynomial residuals")
    bad_refined = bad_mod = 0
    checked = 0
    for n0 in range(4, 15):
        for a0 in range(1, 6):
            R = Q(A, n0, a0, -x)
            V = Q(A, n0, a0 + 1, -x)
            H = Q(A, n0 - 1, a0 + 1, -x)
            S = Q(A, n0 + 1, a0, -x)
            P = Q(A, n0 + 2, a0, -x)
            U = sp.expand(2 * V + x * sp.diff(V, x))
            G0 = sp.expand(
                x * B(V, R) + (n0 + 1) * x * V**2 + V**2 - R * V
            )
            E0 = sp.expand(x * B(S, V) + S * V + (n0 + 1) * a0 * x * V**2)
            rhs0 = sp.expand(a0 * (G0 + n0 * x**2 * B(V, H)))
            bad_refined += int(sp.expand(E0 - rhs0) != 0)
            bad_mod += int(sp.rem(sp.expand(P * V + E0), U, x) != 0)
            checked += 1
    print(
        f"checked={checked}  refined_failures={bad_refined}  "
        f"root_modulo_failures={bad_mod}"
    )

    sep("V6: exact direction sample (m,a)=(6,1)")
    p = Q(A, 6, 1, t)
    s = Q(A, 5, 1, t)
    q1 = Q(A, 5, 2, t)
    v = Q(A, 4, 2, t)
    u = sp.expand(2 * v + t * sp.diff(v, t))
    lam = -sp.Rational(2, 27)
    up = sp.diff(u, t).subs(t, lam)
    wvv = sp.diff(v, t) ** 2 - v * sp.diff(v, t, 2)
    T1 = sp.simplify(-q1.subs(t, lam) * up)
    T2 = sp.simplify(s.subs(t, lam) * up)
    T3 = sp.simplify(
        5 * (2 * v.subs(t, lam) ** 2 + lam**2 * wvv.subs(t, lam))
    )
    total = sp.simplify(T1 + T2 + T3)
    direct = sp.simplify(-p.subs(t, lam) * up)
    print("lambda =", lam)
    print("T1 =", T1, " T2 =", T2, " T3 =", T3)
    print("T1+T2+T3 =", total, " direct -p*u' =", direct)
    print("decomposition residual =", sp.simplify(total - direct))
    print("false shortcut -3v/t =", sp.simplify(-3 * v.subs(t, lam) / lam))
    print("actual u' =", up)

    sep("V7: boundary m=5")
    v5 = Q(A, 3, 2, t)
    u5 = sp.expand(2 * v5 + t * sp.diff(v5, t))
    print("Q_{3,2} =", v5, " UQ_{3,2} =", u5, " roots =", sp.solve(u5, t))

    sep("V8: coupled B(V,H)>0 -- finite evidence only, not the proof")
    bad = 0
    checked = 0
    for n0 in range(4, 15):
        for a0 in range(1, 6):
            V = Q(A, n0, a0 + 1, -x)
            H = Q(A, n0 - 1, a0 + 1, -x)
            w = sp.Poly(B(V, H), x)
            checked += 1
            if w.degree() == 0:
                ok = w.LC() > 0
            else:
                ok = w.LC() > 0 and w.count_roots(-sp.oo, sp.oo) == 0
            bad += int(not ok)
    print(f"checked={checked}  bad={bad} (finite evidence only)")

    sep("V9: first-premise root reduction (7.1) -- exact residual")
    Y0 = 2 * H0 + x * H1
    first_generic = sp.factor(
        sp.rem(
            sp.expand(P0 / a - (a * R0 + x * V1 - n * x * Y0)),
            V0,
            R2,
        )
    )
    print("generic first-premise root residual =", first_generic)
    bad_first = 0
    checked = 0
    for n0 in range(4, 15):
        for a0 in range(1, 6):
            R = Q(A, n0, a0, -x)
            V = Q(A, n0, a0 + 1, -x)
            H = Q(A, n0 - 1, a0 + 1, -x)
            P = Q(A, n0 + 2, a0, -x)
            Y = sp.expand(2 * H + x * sp.diff(H, x))
            res = sp.rem(
                sp.expand(P / a0 - (a0 * R + x * sp.diff(V, x) - n0 * x * Y)),
                V,
                x,
            )
            bad_first += int(res != 0)
            checked += 1
    print(f"actual checked={checked}  root_modulo_failures={bad_first}")

    sep("V10: exact direction sample for the first premise")
    p_first = Q(A, 6, 1, t)
    v_first = Q(A, 4, 2, t)
    lam_first = -sp.Rational(1, 9)
    w_first = sp.simplify(
        -p_first.subs(t, lam_first) * sp.diff(v_first, t).subs(t, lam_first)
    )
    print("lambda =", lam_first, " W(p,v) =", w_first)

    sep("V11: first D2C premise -- finite evidence only, NOT a proof")
    bad = 0
    checked = 0
    for n0 in range(4, 15):
        for a0 in range(1, 6):
            P = Q(A, n0 + 2, a0, -x)
            V = Q(A, n0, a0 + 1, -x)
            w = sp.Poly(B(P, V), x)
            checked += 1
            if w.degree() == 0:
                ok = w.LC() > 0
            else:
                ok = w.LC() > 0 and w.count_roots(-sp.oo, sp.oo) == 0
            bad += int(not ok)
    print(f"checked={checked}  bad={bad} (finite evidence only; all-index status OPEN)")


if __name__ == "__main__":
    main()
