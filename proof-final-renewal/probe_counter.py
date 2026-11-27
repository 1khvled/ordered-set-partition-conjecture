"""Exact counterexamples: renewal is not real-rootedness preserving on the
cone of real-rooted interlacing input families."""
import sympy as sp
import renewal_lib as R

t = R.t


def discrim(p):
    return sp.discriminant(sp.Poly(sp.expand(p), t))


print('=== Exact minimal counterexample, j=3, Laguerre input ===')
P0 = sp.Integer(1)
P1 = 1 + t                      # L_1(-t), root -1
P2 = 1 + 2 * t + t**2 / 2       # L_2(-t), roots -2±sqrt(2)
P = [P0, P1, P2]
Q = [sp.Integer(1)]
for n in [1, 2]:
    rhs, _ = R.renewal_row(P, Q, n, 3)
    Q.append(sp.expand(rhs))
print('P2 =', P2, ' roots:', sorted(float(x) for x in sp.solve(P2, t)))
print('Q2 =', Q[2], ' disc =', discrim(Q[2]))
assert discrim(Q[2]) < 0
print('P1 interlaces P2:', R.interlace([-1], sorted(float(x) for x in sp.solve(P2, t))))

print()
print('=== Degree-correct families: P_n = L_floor(n/(j-1))(-t) ===')
for j in [2, 3, 4]:
    N = 24
    d = lambda n: n // (j - 1)
    P = []
    for n in range(N + 1):
        x = sp.symbols('x')
        p = sp.assoc_laguerre(d(n), 0, x)
        P.append(sp.expand(p.subs(x, -t)))
    Q = [P[0]]
    for n in range(1, N + 1):
        rhs, _ = R.renewal_row(P, Q, n, j)
        Q.append(sp.expand(rhs))
    bad = [n for n in range(N + 1) if not R.real_roots(Q[n])[0]]
    print('j=%d bad n=%s' % (j, bad[:12]))
    for n in bad[:2]:
        print('  Q_%d =' % n, Q[n])
