"""Test candidate all-index coefficient minors between Q_n=O_{n,j} and P_n=O_{n,j-1}."""
import sympy as sp
import renewal_lib as R

t = R.t


def coeff(p, k):
    return sp.Poly(sp.expand(p), t).coeff_monomial(t**k)


for j in [2, 3, 4, 5]:
    Nmax = 20
    P = R.O_rows(j - 1, Nmax)
    Q = R.O_rows(j, Nmax)
    print('=== j', j)
    # Minor between Q_n and P_n: q_{n,k+1} p_{n,k} - q_{n,k} p_{n,k+1}
    bad = []
    for n in range(Nmax + 1):
        dq = sp.Poly(Q[n], t).degree(); dp = sp.Poly(P[n], t).degree()
        for k in range(min(dq, dp)):
            M = coeff(Q[n], k + 1) * coeff(P[n], k) - coeff(Q[n], k) * coeff(P[n], k + 1)
            if M <= 0:
                bad.append((n, k, M))
    print('  Q_n vs P_n minor bad:', bad[:8], 'count', len(bad))

    # Minor Q_n vs P_{n-1}
    bad2 = []
    for n in range(1, Nmax + 1):
        dq = sp.Poly(Q[n], t).degree(); dp = sp.Poly(P[n - 1], t).degree()
        for k in range(min(dq, dp)):
            M = coeff(Q[n], k + 1) * coeff(P[n - 1], k) - coeff(Q[n], k) * coeff(P[n - 1], k + 1)
            if M <= 0:
                bad2.append((n, k))
    print('  Q_n vs P_{n-1} minor bad:', bad2[:8], 'count', len(bad2))

    # Coefficient-wise: q_{n,k} >= p_{n,k}? (Q_n - P_n has (1-t) factor, so no)
    bad3 = [(n, k) for n in range(Nmax + 1)
            for k in range(sp.Poly(Q[n], t).degree() + 1)
            if coeff(Q[n], k) < coeff(P[n], k)]
    print('  q>=p coeffwise bad:', bad3[:8], 'count', len(bad3))
