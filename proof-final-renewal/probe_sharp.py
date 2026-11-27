"""Search for a degree-matched, strictly interlacing input family whose renewal
output fails.  Uses Laguerre rows with alternating parameter."""
import sympy as sp
import renewal_lib as R

x = sp.symbols('x')
t = R.t


def lag(d, alpha, scale=1):
    p = sp.assoc_laguerre(d, alpha, x).subs(x, -t)
    return sp.expand(p * scale)


def roots_sorted(p):
    ok, rr = R.real_roots(p)
    return ok, rr


def try_family(alpha_even, alpha_odd, j=3, nmax=10):
    P = []
    for n in range(nmax + 1):
        d = n // (j - 1)
        alpha = alpha_even if n % 2 == 0 else alpha_odd
        P.append(lag(d, alpha))
    # check consecutive strict interlacing numerically
    ok_all = True
    for n in range(1, nmax + 1):
        ok1, r1 = roots_sorted(P[n - 1])
        ok2, r2 = roots_sorted(P[n])
        if not (ok1 and ok2) or not R.interlace(r1, r2):
            ok_all = False
            break
    Q = [P[0]]
    for n in range(1, nmax + 1):
        rhs, _ = R.renewal_row(P, Q, n, j)
        Q.append(sp.expand(rhs))
    bad = [n for n in range(nmax + 1) if not R.real_roots(Q[n])[0]]
    return ok_all, bad, P, Q


for ae in [0, sp.Rational(1, 2), 1, 2]:
    for ao in [0, sp.Rational(1, 2), 1, 2]:
        if ae == ao:
            continue
        oki, bad, P, Q = try_family(ae, ao, j=3, nmax=10)
        if oki and bad:
            print('alpha_even=%s alpha_odd=%s strict_interlace=%s bad=%s' % (ae, ao, oki, bad[:6]))
            n = bad[0]
            print('  first bad n=%d Q_n=' % n, sp.factor(Q[n]))
