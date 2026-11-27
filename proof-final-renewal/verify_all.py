"""All exact certificates for proof-final-renewal/NOTES.md.

Run with:  py -3.14 verify_all.py
No floating point is used in any certificate; roots are counted with Sturm
sequences and interlacing is certified by the sign of the Wronskian.
"""
import sympy as sp
import renewal_lib as R

t = R.t
x = sp.symbols('x')


def wronskian(f, g):
    return sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))


def exact_interlace(f, g):
    """f,g real-rooted with simple roots; certificate: W(f,g) has no real zero.

    Returns (ok, wronskian_real_roots).
    """
    pf, pg = sp.Poly(f, t), sp.Poly(g, t)
    assert int(sp.count_roots(pf, -sp.oo, sp.oo)) == pf.degree()
    assert int(sp.count_roots(pg, -sp.oo, sp.oo)) == pg.degree()
    W = sp.Poly(wronskian(f, g), t)
    nW = int(sp.count_roots(W, -sp.oo, sp.oo))
    return nW == 0, nW


# ---------------------------------------------------------------- minimal example
print('=== Minimal counterexample j=3, n=2 (unmatched degrees) ===')
P = [sp.Integer(1), 1 + t, 1 + 2 * t + t**2 / 2]
Q = [sp.Integer(1)]
for n in [1, 2]:
    rhs, _ = R.renewal_row(P, Q, n, 3)
    Q.append(sp.expand(rhs))
ok, nW = exact_interlace(P[1], P[2])
disc = int(sp.discriminant(sp.Poly(Q[2], t)))
print('P_1=1+t, P_2=1+2t+t^2/2 ; Wronskian real roots =', nW, 'interlace =', ok)
print('Q_2 =', Q[2], '; discriminant =', disc)
assert ok and disc < 0

# ------------------------------------------------------- degree-matched example
print()
print('=== Degree-matched, strictly interlacing input, j=3, n=8 ===')
j, N = 3, 8


def Pgen(n):
    d = n // (j - 1)
    alpha = sp.Integer(0) if n % 2 == 0 else sp.Rational(1, 2)
    return sp.expand(sp.assoc_laguerre(d, alpha, x).subs(x, -t))


Pm = [Pgen(n) for n in range(N + 1)]
for n in range(N + 1):
    pn = sp.Poly(Pm[n], t)
    assert pn.degree() == n // (j - 1)
    assert int(sp.count_roots(pn, -sp.oo, sp.oo)) == pn.degree()
    assert all(c > 0 for c in pn.all_coeffs())
print('P_n: degree, all-real roots, positive coefficients: OK for n=0..8')
for n in range(1, N + 1):
    ok, nW = exact_interlace(Pm[n - 1], Pm[n])
    assert ok, (n, nW)
print('P_{n-1} interlaces P_n (exact Wronskian certificate): OK for n=1..8')

Qm = [Pm[0]]
for n in range(1, N + 1):
    rhs, _ = R.renewal_row(Pm, Qm, n, j)
    Qm.append(sp.expand(rhs))
q8 = sp.Poly(Qm[8], t)
nreal8 = int(sp.count_roots(q8, -sp.oo, sp.oo))
coeffs8 = q8.all_coeffs()
print('Q_8 =', sp.factor(Qm[8]))
print('  coefficients:', coeffs8)
print('  exact real-root count:', nreal8, '/ degree', q8.degree())
assert any(c < 0 for c in coeffs8)
assert nreal8 < q8.degree()
print('CERTIFIED: renewal output fails real-rootedness.')

# ------------------------------------------------------- renewal identities
print()
print('=== Identity checks ===')
for jj in range(2, 7):
    assert R.verify_renewal(jj, 14)
print('(T4) renewal identity exact for j=2..6, n<=14: OK')
