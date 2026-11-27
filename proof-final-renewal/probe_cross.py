"""Check cross-threshold interlacing and the sign of the convolution at P-roots."""
import sympy as sp
import renewal_lib as R

for j in [2, 3, 4, 5]:
    Nmax = 16
    P = R.O_rows(j - 1, Nmax)
    Q = R.O_rows(j, Nmax)
    print('=== j =', j)
    for n in range(j, Nmax + 1):
        N = n - j + 1
        C = R.binom_conv(Q, P, N)
        okQ, rQ = R.real_roots(Q[n])
        okP, rP = R.real_roots(P[n])
        qvp = R.interlace(rQ, rP) if okQ and okP else None   # Q inside gaps of P
        pvq = R.interlace(rP, rQ) if okQ and okP else None
        # signs of C at roots of P_n
        signs = []
        if okP:
            for mu in rP:
                v = float(C.subs(R.t, mu))
                signs.append(1 if v > 0 else (-1 if v < 0 else 0))
        print('n=%2d dP=%d dQ=%d Q|P=%s P|Q=%s signC@P=%s' %
              (n, sp.Poly(P[n], R.t).degree(), sp.Poly(Q[n], R.t).degree(),
               qvp, pvq, signs))
