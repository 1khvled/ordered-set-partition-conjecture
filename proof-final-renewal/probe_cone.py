"""Test whether the renewal maps an arbitrary real-rooted interlacing family
(P_n) to a real-rooted family (Q_n).  Evidence / counterexample search."""
import sympy as sp
import renewal_lib as R

t = R.t


def laguerre(n, alpha=0):
    """Monic-ish Laguerre L_n^{(alpha)}(-t) has real negative roots."""
    x = sp.symbols('x')
    p = sp.assoc_laguerre(n, alpha, x)
    return sp.expand(p.subs(x, -t))


def trial(Pgen, label, j=2, nmax=12):
    N = nmax + 1
    P = [sp.expand(Pgen(n)) for n in range(N)]
    Q = [P[0]]
    for n in range(1, N):
        rhs, _ = R.renewal_row(P, Q, n, j)
        Q.append(sp.expand(rhs))
    bad = []
    for n in range(N):
        ok, rr = R.real_roots(Q[n])
        if not ok:
            bad.append(n)
    print('%-40s j=%d bad=%s' % (label, j, bad[:8]))
    return bad


# Laguerre family: consecutive Laguerres interlace, all roots negative.
trial(lambda n: laguerre(n, 0), 'Laguerre alpha=0', j=2, nmax=14)
trial(lambda n: laguerre(n, 0), 'Laguerre alpha=0', j=3, nmax=14)
trial(lambda n: laguerre(n, sp.Rational(1, 2)), 'Laguerre alpha=1/2', j=2, nmax=14)
trial(lambda n: laguerre(n, sp.Rational(1, 2)), 'Laguerre alpha=1/2', j=3, nmax=14)

# shifted geometric: P_n = (t+a_n)^n with a_n increasing slowly (interlacing-ish)
for a in [1, 2, sp.Rational(1, 2)]:
    trial(lambda n, a=a: sp.expand((t + a) ** n), 'power (t+%s)^n' % a, j=2, nmax=12)

# Fubini (the actual j=1 family) to confirm it works
trial(lambda n: R.O_rows(n, n)[n] if False else sp.Integer(1), 'dummy', j=2, nmax=3)
