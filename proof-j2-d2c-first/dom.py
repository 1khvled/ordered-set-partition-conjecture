"""Structure of F := x(V')^2 - a B(V,R) modulo V."""
import sympy as sp
from explore import A, Q, B, t

x = sp.Symbol("x", positive=True)


def sep(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)

sep("F = x(V')^2 - aB(V,R); remainder mod V")
for (m, a0) in [(6,1),(7,1),(8,1),(9,1),(7,2),(8,2),(9,2),(10,3)]:
    n = m - 2
    R = Q(A, n, a0, -x)
    V = Q(A, n, a0 + 1, -x)
    Vp = sp.diff(V, x)
    F = sp.expand(x * Vp**2 - a0 * B(V, R))
    rem = sp.rem(F, V, x)
    print(f"\nm={m} a={a0}")
    print("  F =", sp.factor(F))
    print("  rem =", sp.factor(rem))
    rp = sp.Poly(rem, x)
    if rp.degree() > 0:
        print("  rem nreal =", rp.count_roots(-sp.oo, sp.oo), "LC =", rp.LC())
