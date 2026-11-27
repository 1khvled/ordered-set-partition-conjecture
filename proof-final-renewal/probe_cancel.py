"""Check: the renewal cancels the top coefficients of P_n above deg Q_n."""
import sympy as sp
import renewal_lib as R

for j in [3, 4, 5]:
    Nmax = 18
    P = R.O_rows(j - 1, Nmax)
    Q = R.O_rows(j, Nmax)
    print('=== j', j)
    for n in range(j, Nmax + 1):
        N = n - j + 1
        C = R.binom_conv(Q, P, N)
        dP = sp.Poly(P[n], R.t).degree()
        dQ = sp.Poly(Q[n], R.t).degree()
        dC = sp.Poly(C, R.t).degree()
        # coefficients of P_n + (1-t)c_n C above dQ must vanish
        cn = sp.binomial(n, j - 1)
        poly = sp.Poly(P[n] + (1 - R.t) * cn * C, R.t)
        topbad = [k for k in range(dQ + 1, dP + 1) if poly.coeff_monomial(R.t**k) != 0]
        print('n=%2d dP=%d dQ=%d dC=%d cancel_ok=%s' %
              (n, dP, dQ, dC, not topbad))
