"""Find the polynomial G with F = x(V')^2 - aB(V,R) = c x G^2."""
import sympy as sp
from explore import A, Q, B, t

x = sp.Symbol("x", positive=True)


def sep(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)


sep("F = c x G^2 : identify G")
for (m, a0) in [(6,1),(7,1),(8,1),(9,1),(10,1),(11,1),(6,2),(7,2),(8,2),(9,2),(10,3)]:
    n = m - 2
    R = Q(A, n, a0, -x)
    V = Q(A, n, a0 + 1, -x)
    Vp = sp.diff(V, x)
    F = sp.expand(x * Vp**2 - a0 * B(V, R))
    sf = sp.sqf_list(F)
    G = None
    for f, e in sf[1]:
        if e == 2:
            G = f
    c = sp.cancel(F / (x * G**2)) if G is not None else None
    print(f"\nm={m} a={a0}: V={sp.factor(V)}  V'={sp.factor(Vp)}")
    print("   G =", sp.factor(G) if G is not None else "(constant)", " c =", c)
    if G is not None:
        # try to match G against linear combos of V', V/x, R', R/x
        for label, cand in [
            ("V'", Vp),
            ("V/x", sp.cancel(V/x)),
            ("R'", sp.diff(R, x)),
            ("R/x", sp.cancel(R/x)),
            ("a V'+(a+1)R/x", sp.cancel(a0*Vp + (a0+1)*R/x)),
            ("V'+(a+1)R/x", sp.cancel(Vp + (a0+1)*R/x)),
        ]:
            rr = sp.cancel(G / cand) if cand != 0 else None
            if rr is not None and rr.is_number:
                print(f"   *** G = ({rr}) * [{label}]")
