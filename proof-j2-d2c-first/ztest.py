"""Test Z = x^2 R''/R at roots of V and the dominance x(V')^2 > a B(V,R)."""
import sympy as sp
from explore import A, Q, B, t

x = sp.Symbol("x", positive=True)


def sep(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)


sep("Z = x^2 R''/R at roots of V; compare with 3a+2 and 2(a+1)")
for (m, a0) in [(6,1),(7,1),(8,1),(9,1),(10,1),(11,1),
                (6,2),(7,2),(8,2),(9,2),(10,3),(12,1)]:
    n = m - 2
    R = Q(A, n, a0, -x)
    V = Q(A, n, a0 + 1, -x)
    roots = sp.nroots(V, n=40, maxsteps=300)
    print(f"\nm={m} a={a0}: thresholds 3a+2={3*a0+2}, 2(a+1)={2*(a0+1)}")
    for rv in roots:
        if abs(sp.im(rv)) > 1e-18 or sp.re(rv) <= 0:
            continue
        xv = sp.re(rv)
        rv_ = R.subs(x, xv)
        rpp = sp.diff(R, x, 2).subs(x, xv)
        Z = xv**2 * rpp / rv_
        vv = sp.diff(V, x).subs(x, xv)
        bvr = -vv * rv_
        print(f"  x0={float(xv):.6f}  Z={float(Z):.4f}  "
              f"x(V')^2/(aB)={float(xv*vv**2/(a0*bvr)):.4f}  "
              f"dominates? {bool(xv*vv**2 > a0*bvr)}")
