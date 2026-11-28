"""CORRECTED: F = x(V')^2 - a B_x(V,R) with a single x symbol."""
import sympy as sp

x = sp.Symbol("x", positive=True)
t = sp.Symbol("t", real=True)


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
    v = 1
    for k in range(r):
        v *= a + k
    return v


def Q(A, mm, aa, z):
    return sp.expand(sum(rising(aa, r) * A[mm][r] * z**(r-1)
                         for r in range(1, len(A[0])) if A[mm][r]))


def B(f, g, z=x):
    return sp.expand(f * sp.diff(g, z) - sp.diff(f, z) * g)


A = astirling_table(26)


def sep(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)

sep("CORRECTED F = x(V')^2 - aB(V,R), and Delta at roots")
for (m, a0) in [(6,1),(7,1),(8,1),(9,1),(10,1),(11,1),(7,2),(8,2),(9,2),(10,3)]:
    n = m - 2
    R = Q(A, n, a0, -x)
    V = Q(A, n, a0 + 1, -x)
    Vp = sp.diff(V, x)
    H = Q(A, n - 1, a0 + 1, -x)
    Y = sp.expand(2*H + x*sp.diff(H, x))
    F = sp.expand(x * Vp**2 - a0 * B(V, R))
    Delta = sp.expand(x*Vp**2 + n*x*B(V,Y) - a0*B(V,R))
    # evaluate F and Delta at roots
    roots = sp.nroots(V, n=30, maxsteps=200)
    vals = []
    for rv in roots:
        if abs(sp.im(rv)) < 1e-18 and sp.re(rv) > 0:
            xv = sp.re(rv)
            fv = F.subs(x, xv)
            dv = Delta.subs(x, xv)
            bvr = B(V,R).subs(x,xv)
            vals.append((float(xv), float(fv), float(bvr), float(dv)))
    print(f"\nm={m} a={a0}: n={n}")
    print("  F factor =", sp.factor(F))
    print("  Delta factor =", sp.factor(Delta))
    for xv, fv, bvr, dv in vals:
        print(f"    x0={xv:.5f}: F={fv:.4f} aB={a0*bvr:.4f} Delta={dv:.4f}")
