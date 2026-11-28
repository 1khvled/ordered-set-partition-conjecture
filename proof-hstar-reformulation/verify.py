#!/usr/bin/env python3
"""Exact audits for the v/h* reformulation (no builds, no floating point).

Every algebraic residual used in NOTES.md is printed.  Finite ranges are
consistency checks; the all-index arguments are written separately in NOTES.md.
"""
from __future__ import annotations

import math
from fractions import Fraction

import sympy as sp

v, y, z, m, t, x = sp.symbols("v y z m t x")
JMAX = 8
NMAX = 12


def mul_series(a: list[Fraction], b: list[Fraction], nmax: int) -> list[Fraction]:
    out = [Fraction(0) for _ in range(nmax + 1)]
    for i, ai in enumerate(a):
        if not ai:
            continue
        for k, bk in enumerate(b[: nmax + 1 - i]):
            if bk:
                out[i + k] += ai * bk
    return out


def inv_series(a: list[Fraction], nmax: int) -> list[Fraction]:
    assert a[0] == 1
    out = [Fraction(0) for _ in range(nmax + 1)]
    out[0] = Fraction(1)
    for n in range(1, nmax + 1):
        out[n] = -sum(a[k] * out[n - k] for k in range(1, n + 1))
    return out


def ordered_polys(j: int, nmax: int) -> list[sp.Expr]:
    """Original first-block recurrence, independent of the universal one."""
    rows = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        p = sp.Integer(0)
        for r in range(1, min(n, j - 1) + 1):
            p += math.comb(n, r) * rows[n - r]
        for r in range(j, n + 1):
            p += math.comb(n, r) * t * rows[n - r]
        rows.append(sp.expand(p))
    return rows


def definition_polys(j: int, nmax: int) -> list[sp.Expr]:
    """Invert 1-A-tC directly as a truncated formal EGF series."""
    den = [sp.Integer(0) for _ in range(nmax + 1)]
    den[0] = sp.Integer(1)
    for r in range(1, nmax + 1):
        den[r] = -sp.Rational(1, math.factorial(r)) * (t if r >= j else 1)
    f = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        f.append(sp.expand(-sum(den[k] * f[n - k] for k in range(1, n + 1))))
    return [sp.expand(math.factorial(n) * f[n]) for n in range(nmax + 1)]


def hstar_from_egf(j: int, n: int, mv: int) -> sp.Integer:
    """n![z^n] (1-A)^mv/(2-e^z)^(mv+1), directly as a truncated series."""
    nmax = n
    a = [Fraction(0) for _ in range(nmax + 1)]
    a[0] = Fraction(1)
    for r in range(1, min(j, nmax + 1)):
        a[r] -= Fraction(1, math.factorial(r))
    e = [Fraction(0) for _ in range(nmax + 1)]
    e[0] = Fraction(1)
    for r in range(1, nmax + 1):
        e[r] = -Fraction(1, math.factorial(r))
    einv = inv_series(e, nmax)
    deninv = [Fraction(1)] + [Fraction(0)] * nmax
    for _ in range(mv + 1):
        deninv = mul_series(deninv, einv, nmax)
    numerator = [Fraction(1)] + [Fraction(0)] * nmax
    for _ in range(mv):
        numerator = mul_series(numerator, a, nmax)
    coeff = mul_series(numerator, deninv, nmax)[n]
    value = coeff * math.factorial(n)
    assert value.denominator == 1
    return sp.Integer(value.numerator)


def c_series_from_egf(j: int, nmax: int) -> list[list[Fraction]]:
    """c[n][r] = n![z^n] C^r/E^(r+1), with C=e^z-1-A."""
    e = [Fraction(0) for _ in range(nmax + 1)]
    e[0] = Fraction(1)
    for q in range(1, nmax + 1):
        e[q] = -Fraction(1, math.factorial(q))
    c = [Fraction(0) for _ in range(nmax + 1)]
    for q in range(j, nmax + 1):
        c[q] = Fraction(1, math.factorial(q))
    einv = inv_series(e, nmax)
    c_over_e = mul_series(c, einv, nmax)
    hpow = [Fraction(1)] + [Fraction(0)] * nmax
    out = [[Fraction(0) for _ in range(nmax // j + 1)] for _ in range(nmax + 1)]
    for r in range(nmax // j + 1):
        prod = mul_series(einv, hpow, nmax)
        for n in range(nmax + 1):
            out[n][r] = prod[n] * math.factorial(n)
        hpow = mul_series(hpow, c_over_e, nmax)
    return out


def binom_poly(arg: sp.Expr, r: int) -> sp.Expr:
    return sp.prod((arg - q for q in range(r)), sp.Integer(1)) / sp.Integer(math.factorial(r))


def finite_difference_coeffs(g: sp.Expr, degree: int) -> list[sp.Expr]:
    return [
        sp.simplify(
            sum((-1) ** (r - q) * math.comb(r, q) * g.subs(m, q) for q in range(r + 1))
        )
        for r in range(degree + 1)
    ]


def gamma_coeffs(g: sp.Expr, degree: int) -> list[sp.Expr]:
    aa = sp.symbols(f"a0:{degree + 1}")
    form = sum(aa[k] * m**k * (1 - m) ** (degree - k) for k in range(degree + 1))
    equations = sp.Poly(sp.expand(form - g), m).coeffs()
    solution = sp.solve(equations, aa, dict=True)[0]
    return [sp.simplify(solution[a]) for a in aa]


def main() -> None:
    checks = 0

    def zero(label: str, residual: sp.Expr) -> None:
        nonlocal checks
        residual = sp.simplify(residual)
        print(f"RESIDUAL {label} = {residual}")
        assert residual == 0, (label, residual)
        checks += 1

    # Generic chain-rule identities used to transform the universal recurrence.
    f = v**3 + 2 * v + 1
    zero(
        "chain_T",
        (t * sp.diff(f.subs(v, t / (1 - t)), t)).subs(t, v / (1 + v))
        - v * (1 + v) * sp.diff(f, v),
    )
    zero(
        "chain_t(t-1)d",
        (t * (t - 1) * sp.diff(f.subs(v, t / (1 - t)), t)).subs(t, v / (1 + v))
        + v * sp.diff(f, v),
    )
    zero("map_v(t)_substituted", (v - t / (1 - t)).subs(t, v / (1 + v)))
    zero("map_y(v)_substituted", (y + 1 / v).subs(v, -1 / y))
    zero("map_y(t)_substituted", (y - (1 - 1 / t)).subs(t, 1 / (1 - y)))
    zero("map_x(t)_substituted", (x - (t - 1)).subs(t, 1 + x))
    zero("map_x(y)_substituted", (x - y / (1 - y)).subs(x, y / (1 - y)))

    # Direction and monotonicity audit, with explicit samples and endpoints.
    print("DIRECTION y(t)=1-1/t, dy/dt=1/t^2 > 0 on t<0 and on 0<t<1")
    for tv in [sp.Rational(-2), sp.Rational(-1), sp.Rational(-1, 2), sp.Rational(1, 2)]:
        yv = sp.simplify(1 - 1 / tv)
        print(f"SAMPLE t={tv} -> y={yv}; y>1 is {yv > 1}; derivative={1 / tv**2}")
    print("ENDPOINTS t->-oo: y->1+; t->0-: y->+oo; t->0+: y->-oo; t->1-: y->0-")
    print("DIRECTION v(t)=t/(1-t), dv/dt=1/(1-t)^2 > 0 for t<0")
    for tv in [sp.Rational(-2), sp.Rational(-1), sp.Rational(-1, 2)]:
        print(f"SAMPLE t={tv} -> v={sp.simplify(tv / (1 - tv))}")
    print("ENDPOINTS t->-oo: v->-1+; t->0-: v->0-")
    print("DIRECTION x(y)=y/(1-y), dx/dy=1/(1-y)^2 > 0 for y>1")
    for yv in [sp.Rational(3, 2), sp.Integer(2), sp.Integer(3)]:
        print(f"SAMPLE y={yv} -> x={sp.simplify(yv / (1 - yv))}")
    print("ENDPOINTS y->1+: x->-oo; y->+oo: x->-1-")

    # Exact generating-function algebra.
    a = sum(z**r / math.factorial(r) for r in range(1, JMAX))
    c = sp.exp(z) - sum(z**r / math.factorial(r) for r in range(JMAX))
    nfac = 1 - a
    efac = 2 - sp.exp(z)
    zero("denominator_t", 1 - a - t * c - ((1 + t) - (1 - t) * a - t * sp.exp(z)))
    zero(
        "denominator_v",
        sp.cancel((1 + v) / (1 - a - v * (sp.exp(z) - 2)) - 1 / (1 - a - t * c)).subs(
            t, v / (1 + v)
        ),
    )
    zero("C=N-E", c - (nfac - efac))
    w = z ** (JMAX - 1) / math.factorial(JMAX - 1)
    dfac = efac - x * c
    hfac = 1 / dfac
    zero("stable-denominator derivative", sp.diff(dfac, z) - (dfac - 2 - x * w))
    zero(
        "stable-denominator Riccati",
        sp.diff(hfac, z) - ((2 + x * w) * hfac**2 - hfac),
    )
    zero(
        "hstar_geometric",
        (1 - y) / (efac - y * nfac)
        - (1 - y) / efac / (1 - y * nfac / efac),
    )
    zero("hstar_at_t", sp.cancel((1 - y) / (efac - y * nfac) - 1 / (1 - a - t * c)).subs(t, 1 / (1 - y)))

    gammabad = []
    crootbad = []
    for j in range(2, JMAX + 1):
        ctable = c_series_from_egf(j, NMAX)
        rows = ordered_polys(j, NMAX)
        defrows = definition_polys(j, NMAX)
        rrows: list[sp.Expr] = []
        wrows: list[sp.Expr] = []
        for n, p in enumerate(rows):
            zero(f"definition j={j} n={n}", p - defrows[n])
            d = n // j
            assert sp.degree(p, t) == d
            assert p.subs(t, 0) > 0
            rn = sp.Poly(sp.cancel((1 + v) ** d * p.subs(t, v / (1 + v))), v).as_expr()
            wn = sp.Poly(sp.cancel((1 - y) ** d * p.subs(t, 1 / (1 - y))), y).as_expr()
            zero(f"R-W j={j} n={n}", (-1) ** d * y**d * rn.subs(v, -1 / y) - wn)
            zero(
                f"W-at-one cancellation j={j} n={n}",
                wn.subs(y, 1) - p.coeff(t, d),
            )
            zero(
                f"R-at-minus-one cancellation j={j} n={n}",
                rn.subs(v, -1) - (-1) ** d * p.coeff(t, d),
            )
            rrows.append(rn)
            wrows.append(wn)

            # g reconstructed from W, then checked directly against its EGF.
            g = sp.expand(
                sum(
                    wn.coeff(y, q) * binom_poly(m + d - q, d)
                    for q in range(d + 1)
                )
            )
            for mv in range(0, d + 5):
                zero(
                    f"g-EGF j={j} n={n} m={mv}",
                    g.subs(m, mv) - hstar_from_egf(j, n, mv),
                )

            # c is the binomial transform of g and has the stable-denominator EGF.
            cs = finite_difference_coeffs(g, d)
            for r in range(d + 1):
                zero(f"c-EGF j={j} n={n} r={r}", cs[r] - sp.Integer(ctable[n][r]))
            zero(
                f"W-c j={j} n={n}",
                wn - sum(cs[r] * y**r * (1 - y) ** (d - r) for r in range(d + 1)),
            )
            cx = sum(cs[r] * x**r for r in range(d + 1))
            zero(
                f"W-c-substitution j={j} n={n}",
                wn - (1 - y) ** d * cx.subs(x, y / (1 - y)),
            )
            zero(
                f"c-shift j={j} n={n}",
                cx - p.subs(t, 1 + x),
            )

            # Finite evidence only: actual g rows happen to be gamma-positive and
            # real-rooted below -1 in this box.  This is not used as a proof.
            if n >= j:
                roots = sp.Poly(g, m).all_roots()
                if not all(r.is_real and r < -1 for r in roots):
                    gammabad.append((j, n, "g-root", roots))
                gc = gamma_coeffs(g, d)
                if not all(a >= 0 for a in gc):
                    gammabad.append((j, n, "gamma", gc))
                croots = sp.Poly(cx, x).all_roots()
                if not all(r.is_real and r < -1 for r in croots):
                    crootbad.append((j, n, croots))

        # Universal recurrence transformed exactly to R.
        for n in range(1, NMAX + 1):
            d = n // j
            dprev = (n - 1) // j
            s = 1 + v if n % j == 0 else 1
            rhs = s * (
                n * rrows[n - 1]
                + 2 * v * (1 + v) * sp.diff(rrows[n - 1], v)
                - 2 * dprev * v * rrows[n - 1]
            )
            for r in range(2, j):
                if n >= r:
                    sr = 1 + v if (n - r) // j != d else 1
                    rhs += math.comb(n, r) * sr * rrows[n - r]
            if n >= j:
                rhs += (
                    math.comb(n - 1, j - 1)
                    * v
                    * (d * rrows[n - j] - (1 + v) * sp.diff(rrows[n - j], v))
                )
            zero(f"R-recurrence j={j} n={n}", rrows[n] - rhs)
        print(f"EXACT j={j}: n=0..{NMAX}, definition, universal/R recurrence, W/c/g maps")

    print(f"FINITE g gamma/root box: failures={len(gammabad)}; details={gammabad}")
    print(f"FINITE c root box: failures={len(crootbad)}; details={crootbad}")

    # Positive-root trap survives in v: j=2,n=5 lag term.
    rows2 = ordered_polys(2, 5)
    q = rows2[3]
    tq = sp.expand(t * (q + (t - 1) * sp.diff(q, t)))
    normalized = sp.Poly(
        sp.cancel((1 + v) ** 2 * tq.subs(t, v / (1 + v))), v
    ).as_expr()
    zero("trap normalization", normalized - (13 * v**2 - v))
    print(f"TRAP O_3={q}; tQ(O_3)={tq}; roots={sp.solve(tq, t)}")
    print(f"TRAP normalized={normalized}; roots={sp.solve(normalized, v)}")
    assert sp.solve(tq, t) == [0, sp.Rational(1, 14)]
    assert sp.solve(normalized, v) == [0, sp.Rational(1, 13)]

    # Exact obstruction supplied in the task, audited in both bases.
    gmodel = (m + 2) ** 2
    ccoef = finite_difference_coeffs(gmodel, 2)
    wmodel = sp.expand(sum(ccoef[r] * y**r * (1 - y) ** (2 - r) for r in range(3)))
    gen_model = sp.cancel(
        y * (1 + y) / (1 - y) ** 3
        + 4 * y / (1 - y) ** 2
        + 4 / (1 - y)
    )
    zero("obstruction generating function", (1 - y) ** 3 * gen_model - wmodel)
    zero(
        "obstruction g-binomial transform",
        gmodel - sum(ccoef[r] * binom_poly(m, r) for r in range(3)),
    )
    gammamodel = gamma_coeffs(gmodel, 2)
    cmodel = sum(ccoef[r] * x**r for r in range(3))
    zero("obstruction gamma transform", gmodel - sum(gammamodel[k] * m**k * (1 - m) ** (2 - k) for k in range(3)))
    print(f"OBSTRUCTION g(m)={gmodel}; roots={sp.solve(gmodel, m)}; gamma={gammamodel}")
    print(f"OBSTRUCTION c(x)={cmodel}; c discriminant={sp.discriminant(cmodel, x)}")
    print(f"OBSTRUCTION W(y)={wmodel}; W discriminant={sp.discriminant(wmodel, y)}")
    assert gammamodel == [4, 12, 9]
    assert ccoef == [4, 5, 2]
    assert sp.discriminant(cmodel, x) == -7
    assert sp.discriminant(wmodel, y) == -7

    # Boundary/degeneracy enumeration.
    print("BOUNDARY n=0: d=0, O=1, R=W=c=1; no nonconstant roots")
    for j in range(2, JMAX + 1):
        print(
            f"BOUNDARY j={j}: n<j gives d=0; n=j gives d=1; "
            f"n=2j gives d=2 and the (1+v) normalization jump"
        )
    print("BOUNDARY inactive r: c[n][r]=0 for r>floor(n/j); g degree is exactly d")
    print("BOUNDARY W(1)=leading coefficient of O > 0; R(-1)=(-1)^d times it; neither point is a root")
    print("BOUNDARY y=0 means t=1; O_n(1)>0, so y=0 is not a root")
    print(f"TOTAL exact identity checks = {checks}")


if __name__ == "__main__":
    main()
