"""Exact residuals and direction checks for the first D2C premise study.

No installs/builds. SymPy (already installed) only. Finite checks labelled.
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
    v = 1
    for k in range(r):
        v *= a + k
    return v


def Q(A, m, a, z):
    return sp.expand(sum(rising(a, r) * A[m][r] * z**(r-1)
                         for r in range(1, len(A[0])) if A[m][r]))


def W(f, g, z=t):
    return sp.expand(sp.diff(f, z) * g - f * sp.diff(g, z))


def B(f, g, z=x):
    return sp.expand(f * sp.diff(g, z) - sp.diff(f, z) * g)


def sep(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)


A = astirling_table(24)

sep("V1: (M), (V), (R), (7.1) residuals")
bad = 0
checked = 0
for m in range(4, 16):
    for a0 in range(1, 6):
        q = Q(A, m, a0, t)
        qm1 = Q(A, m - 1, a0, t)
        qm2 = Q(A, m - 2, a0 + 1, t)
        qa1 = Q(A, m, a0 + 1, t)
        bad += int(sp.expand(q - (sp.diff(t*qm1, t) + (m-1)*a0*t*qm2)) != 0)
        bad += int(sp.expand(qa1 - q - sp.diff(t*q, t)/a0) != 0)
        checked += 1
print(f"checked={checked} bad={bad}")

sep("V2: coefficient-ratio monotonicity (proof of coefficient positivity)")
bad = 0
total = 0
for m in range(6, 21):
    for a0 in range(1, 7):
        for r in range(2, m // 2 + 1):
            for s in range(1, r):
                if A[m][r] == 0 or A[m-2][s] == 0:
                    continue
                lhs = (a0+s) * A[m][r] * A[m-2][s]
                rhs = (a0+r) * A[m][s] * A[m-2][r]
                total += 1
                bad += int(lhs < rhs)
print(f"pairs={total} violations={bad}")

sep("V3: all coefficients of W(p,w) positive (finite evidence)")
bad = 0
checked = 0
for m in range(6, 20):
    for a0 in range(1, 7):
        p = Q(A, m, a0, t)
        w = Q(A, m - 2, a0 + 1, t)
        co = sp.Poly(W(p, w), t).all_coeffs()
        bad += int(not all(c > 0 for c in co))
        checked += 1
print(f"checked={checked} bad={bad}")

sep("V4: exact identity (7.1) mod V, and exact (7.3) reduction")
# generic residual of P/a = aR + xV' - nxY modulo V
a, n = sp.symbols("a n", positive=True)
R0, R1, R2 = sp.symbols("R0 R1 R2")
H0, H1 = sp.symbols("H0 H1")
V0 = ((a+1)*R0 + x*R1)/a
V1 = ((a+2)*R1 + x*R2)/a
Y0 = 2*H0 + x*H1
S0 = sp.expand(a*(V0 - R0 - n*x*H0))
S1 = sp.expand(a*(V1 - R1 - n*(H0 + x*H1)))
P0 = sp.expand(S0 + x*S1 - (n+1)*a*x*V0)
res71 = sp.rem(sp.expand(P0/a - (a*R0 + x*V1 - n*x*Y0)), V0, R2)
print("generic (7.1) residual =", sp.factor(res71))

bad = 0
checked = 0
for m in range(6, 16):
    for a0 in range(1, 6):
        nn = m - 2
        R = Q(A, nn, a0, -x)
        V = Q(A, nn, a0 + 1, -x)
        H = Q(A, nn - 1, a0 + 1, -x)
        Y = sp.expand(2*H + x*sp.diff(H, x))
        P = Q(A, nn + 2, a0, -x)
        res = sp.rem(sp.expand(P/a0 - (a0*R + x*sp.diff(V,x) - nn*x*Y)), V, x)
        bad += int(res != 0)
        checked += 1
print(f"actual (7.1): checked={checked} bad={bad}")

sep("V5: dominance criterion Z < 3a+2 (exact equivalence)")
print("Needs exact verification; see NOTES §4. Finite evidence below:")
for (m, a0) in [(6,1),(8,1),(10,1),(8,2),(10,3)]:
    nn = m - 2
    R = Q(A, nn, a0, -x)
    V = Q(A, nn, a0 + 1, -x)
    rv = [r for r in sp.nroots(V, n=30, maxsteps=200)
          if abs(sp.im(r)) < 1e-18 and sp.re(r) > 0]
    print(f"  m={m} a={a0}: 3a+2={3*a0+2}, Z at roots = "
          + str([round(float(sp.re(r)**2*sp.diff(R,x,2).subs(x,sp.re(r))/R.subs(x,sp.re(r))), 3)
                 for r in rv]))

sep("V6: direction sample (m,a)=(6,1) -- first premise")
p = Q(A, 6, 1, t)
v = Q(A, 4, 2, t)
lam = -sp.Rational(1, 9)
wv = sp.simplify(-p.subs(t, lam)*sp.diff(v, t).subs(t, lam))
print("lambda =", lam, " W(p,v) =", wv)
