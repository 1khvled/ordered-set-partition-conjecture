"""Probe: does the bilinear convolution (Q star P) have real negative roots,
and does it interlace P_n / Q_n?  Evidence only."""
import sympy as sp
import renewal_lib as R

for j in [2, 3, 4, 5]:
    Nmax = 18
    P = R.O_rows(j - 1, Nmax)
    Q = R.O_rows(j, Nmax)
    print('=== j =', j)
    for n in range(j, Nmax + 1):
        N = n - j + 1
        C = R.binom_conv(Q, P, N)
        okC, rC = R.real_roots(C)
        okP, rP = R.real_roots(P[n])
        okQ, rQ = R.real_roots(Q[n])
        # check signs at 0
        def sgn(p):
            v = sp.expand(p).subs(R.t, 0)
            return sp.sign(v)
        negC = all(x <= 1e-9 for x in rC)
        interCP = R.interlace(rC, rP) if okC and okP else None
        interCQ = R.interlace(rC, rQ) if okC and okQ else None
        interPQ = R.interlace(rP, rQ) if okP and okQ else None
        print('n=%2d degQ=%2d Crr=%s negC=%s C|P=%s C|Q=%s P|Q=%s' %
              (n, sp.Poly(Q[n], R.t).degree(), okC, negC, interCP, interCQ, interPQ))
