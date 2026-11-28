"""Probe the first-premise inequality at exact roots of V."""
import sympy as sp
from explore import A, Q, W, B, t, x

x = sp.symbols("x", positive=True)


def sep(s):
    print("\n" + "=" * 68)
    print(s)
    print("=" * 68)


sep("A. Positivity of W(Q_{m,a}, Q_{m-2,a}) (same a, distance 2)")
for a0 in [1, 2]:
    for m in range(6, 12):
        p = Q(A, m, a0, t)
        h = Q(A, m - 2, a0, t)
        w = W(p, h)
        pol = sp.Poly(w, t)
        nreal = pol.count_roots(-sp.oo, sp.oo) if pol.degree() > 0 else 0
        print(f"m={m} a={a0}: deg={pol.degree()} LC={pol.LC()} nreal={nreal} "
              f"const={pol.eval(0)}")

sep("B. Root-of-V data for (7.3): P V'/a decomposition")
for (m, a0) in [(6, 1), (7, 1), (8, 1), (7, 2), (8, 2), (9, 1)]:
    n = m - 2
    R = Q(A, n, a0, -x)
    V = Q(A, n, a0 + 1, -x)
    H = Q(A, n - 1, a0 + 1, -x)
    Y = sp.expand(2 * H + x * sp.diff(H, x))
    P = Q(A, n + 2, a0, -x)
    roots = sp.nroots(V, n=30, maxsteps=200)
    print(f"\nm={m} a={a0}, V={sp.factor(V)}")
    print("  roots of V (x>0):", [str(r) for r in roots if abs(sp.im(r)) < 1e-20 and sp.re(r) > 0])
    for rv in roots:
        if abs(sp.im(rv)) > 1e-20 or sp.re(rv) <= 0:
            continue
        xv = sp.re(rv)
        vv = sp.diff(V, x).subs(x, xv)
        rv_ = R.subs(x, xv)
        yv_ = Y.subs(x, xv)
        bvr = B(V, R).subs(x, xv)
        bvy = B(V, Y).subs(x, xv)
        pv = P.subs(x, xv)
        delta = xv * vv**2 + n * xv * bvy - a0 * bvr
        print(f"  x0={float(xv):.6f}: V'={float(vv):.4f} R={float(rv_):.4f} "
              f"Y={float(yv_):.4f} B(V,R)={float(bvr):.4f} B(V,Y)={float(bvy):.4f} "
              f"Delta={float(delta):.6f}  PV'/a={float(pv*vv/a0):.6f}")

sep("C. Coefficient positivity of W(p,w)")
for a0 in [1, 2, 3]:
    for m in range(6, 14):
        p = Q(A, m, a0, t)
        w = Q(A, m - 2, a0 + 1, t)
        co = sp.Poly(W(p, w), t).all_coeffs()
        print(f"m={m} a={a0}: all coeffs >0? {all(c>0 for c in co)}")
