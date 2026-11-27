"""Exact certificates for the sharp counterexample in NOTES.md.

Input family (j=3): P_n(t) = L_{floor(n/2)}^{(alpha_n)}(-t),
alpha_n = 0 (n even), 1/2 (n odd),  n = 0..8.
Renewal output Q_8 has a negative t-coefficient, hence is not real-rooted
with nonpositive roots.
"""
import sympy as sp
import renewal_lib as R

t = R.t
x = sp.symbols('x')
j = 3
N = 8


def Pgen(n):
    d = n // (j - 1)
    alpha = sp.Integer(0) if n % 2 == 0 else sp.Rational(1, 2)
    return sp.expand(sp.assoc_laguerre(d, alpha, x).subs(x, -t))


P = [Pgen(n) for n in range(N + 1)]
Q = [P[0]]
for n in range(1, N + 1):
    rhs, _ = R.renewal_row(P, Q, n, j)
    Q.append(sp.expand(rhs))

print('--- P rows ---')
for n in range(N + 1):
    pp = sp.Poly(P[n], t)
    nreal = int(sp.count_roots(pp, -sp.oo, sp.oo))
    # positive coefficients => no positive roots
    pos = all(c > 0 for c in pp.all_coeffs())
    print('P_%d deg=%d real=%d poscoeff=%s' % (n, pp.degree(), nreal, pos))
    assert nreal == pp.degree() and pos

print('--- consecutive strict interlacing (Wronskian sign / root alternation) ---')
for n in range(1, N + 1):
    a, b = P[n - 1], P[n]
    # exact: roots of the two polynomials strictly alternate iff the Wronskian
    # has constant sign and degrees differ by 0 or 1 (plus leading sign check)
    ok, ra = R.real_roots(a)
    okb, rb = R.real_roots(b)
    inter = R.interlace(ra, rb)
    print('P_%d | P_%d : %s  roots %s / %s' % (n - 1, n, inter,
          [round(v, 4) for v in ra], [round(v, 4) for v in rb]))
    assert inter

print('--- renewal output Q_8 ---')
q8 = sp.Poly(Q[8], t)
print('Q_8 =', sp.factor(Q[8]))
coeffs = q8.all_coeffs()
print('coeffs (high->low):', coeffs)
nreal8 = int(sp.count_roots(q8, -sp.oo, sp.oo))
print('degree =', q8.degree(), ' exact # real roots =', nreal8)
assert any(c < 0 for c in coeffs), 'expected a negative coefficient'
assert nreal8 < q8.degree(), 'expected nonreal roots'
print('CERTIFIED: Q_8 has negative coefficients and only', nreal8, 'real roots of', q8.degree())
