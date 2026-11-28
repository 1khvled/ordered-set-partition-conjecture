"""Test: real-rooted negative-zero polys with positive coeffs and monotone
coefficient ratios -- does that force W(f,g)>0 (interlacing)?"""
import random
import sympy as sp

t = sp.Symbol("t", real=True)

random.seed(12345)


def poly_from_roots(rs):
    p = sp.Integer(1)
    for r in rs:
        p = sp.expand(p * (t + r))
    return p


def ratio_increasing(f, g):
    pf = sp.Poly(f, t)
    pg = sp.Poly(g, t)
    df, dg = pf.degree(), pg.degree()
    N = min(df, dg)
    vals = []
    for k in range(N + 1):
        a = pf.coeff_monomial(t ** k)
        b = pg.coeff_monomial(t ** k)
        vals.append(sp.Rational(a, b))
    return all(vals[i] < vals[i + 1] for i in range(N)), vals


def W(f, g):
    return sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))


found = 0
tested = 0
for trial in range(20000):
    nf = random.randint(2, 4)
    ng = nf - 1
    rs_f = [sp.Rational(random.randint(1, 40), random.randint(1, 10)) for _ in range(nf)]
    rs_g = [sp.Rational(random.randint(1, 40), random.randint(1, 10)) for _ in range(ng)]
    f = poly_from_roots(rs_f)
    g = poly_from_roots(rs_g)
    mono, vals = ratio_increasing(f, g)
    if not mono:
        continue
    tested += 1
    w = W(f, g)
    pol = sp.Poly(w, t)
    ok = pol.LC() > 0 and (pol.degree() == 0 or pol.count_roots(-sp.oo, sp.oo) == 0)
    if not ok:
        found += 1
        if found <= 5:
            print("COUNTEREXAMPLE")
            print("  f =", f, " roots", rs_f)
            print("  g =", g, " roots", rs_g)
            print("  ratios", [float(v) for v in vals])
            print("  W =", sp.factor(w))
print(f"tested monotone ratio pairs: {tested}, failures: {found}")
