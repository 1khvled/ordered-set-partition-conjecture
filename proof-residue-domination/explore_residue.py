"""Exact exploration of the residue-mass domination G_q(mu)>0.

Builds the ordered big-block rows O_n(t) over Q[t] by the first-block
recurrence, the inner rows H_n(t)=O_n - sum_{r=1}^{j-1} C(n,r) O_{n-r},
and tests, for q=O_{n-d}:

  (I)  exact polynomial identity  W(H_n,q) + t q' S = 0  mod q,
       where S = sum_{r=j}^n C(n,r) O_{n-r}  (so G_q = S/(t q') at a root).
  (II) exact identity  H_n = t S  (polynomial equality).
  (III) residue values G_q(mu) and aligned interlacing of H_n with q.
  (IV) far-lag residue sign variation (termwise positivity fails).

No floats are used as certificates; (I)/(II) are exact symbolic, and
root-based statements use sympy real_roots (algebraic isolation) where
feasible, else high-precision nroots with an explicit warning.
"""
import sympy as sp

t = sp.symbols('t')


def rows(j, nmax):
    """O_n and H_n for n=0..nmax, exact in Q[t]."""
    O = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = sp.Integer(0)
        for r in range(1, min(j - 1, n) + 1):
            s += sp.binomial(n, r) * O[n - r]
        for r in range(j, n + 1):
            s += t * sp.binomial(n, r) * O[n - r]
        O.append(sp.expand(s))
    H = []
    for n in range(nmax + 1):
        s = O[n]
        for r in range(1, min(j - 1, n) + 1):
            s -= sp.binomial(n, r) * O[n - r]
        H.append(sp.expand(s))
    return O, H


def W(p, q):
    return sp.expand(sp.diff(p, t) * q - p * sp.diff(q, t))


def check_identities(jmax=5, nmax=22):
    """Exact checks (I),(II) over a moderate grid."""
    tested = 0
    for j in range(2, jmax + 1):
        O, H = rows(j, nmax)
        for n in range(j, nmax + 1):
            S = sum(sp.binomial(n, r) * O[n - r] for r in range(j, n + 1))
            S = sp.expand(S)
            assert sp.expand(H[n] - t * S) == 0, ('II', j, n)
            for d in range(1, j + 1):
                q = O[n - d]
                if sp.degree(q, t) < 1:
                    continue
                # W(H_n,q)(mu) = -mu q'(mu) S(mu) at roots mu of q
                poly = sp.expand(W(H[n], q) + t * sp.diff(q, t) * S)
                rem = sp.rem(sp.Poly(poly, t), sp.Poly(q, t))
                assert rem.as_expr() == 0, ('I', j, n, d)
                tested += 1
    print('identity checks passed:', tested)


def residue_and_interlacing(j, n, d):
    """High-precision residue values and aligned interlacing for one triple."""
    O, H = rows(j, n)
    q = O[n - d]
    if sp.degree(q, t) < 1:
        return None
    S = sum(sp.binomial(n, r) * O[n - r] for r in range(j, n + 1))
    # G_q(mu) = S(mu)/q'(mu);  H_n=tS gives H_n(mu)=mu*S(mu), so also
    # G_q(mu) = H_n(mu)/(mu q'(mu)).  Must NOT divide S by mu here.
    qroots = [sp.re(r) for r in sp.nroots(q, n=40) if abs(sp.im(r)) < 1e-25]
    qroots = sorted([r for r in qroots if r < 0])
    if not qroots:
        return None
    Gvals = []
    for mu in qroots:
        qp = sp.diff(q, t).subs(t, mu)
        Sv = S.subs(t, mu)
        Gvals.append(sp.N(Sv / qp, 20))
    # roots of H_n
    Hroots = sorted([sp.re(r) for r in sp.nroots(H[n], n=40)
                     if abs(sp.im(r)) < 1e-25 and sp.re(r) < 0])
    return qroots, Gvals, Hroots


def report():
    print('=== residue values / interlacing samples ===')
    for (j, n, d) in [(2, 8, 2), (2, 10, 2), (3, 12, 3), (3, 9, 3),
                      (4, 16, 4), (2, 9, 1), (3, 11, 2)]:
        r = residue_and_interlacing(j, n, d)
        if r is None:
            continue
        qroots, Gvals, Hroots = r
        print(f'j={j} n={n} d={d}: deg q={len(qroots)}, deg H={len(Hroots)}')
        print('  G values:', [float(g) for g in Gvals])
        print('  q roots :', [float(x) for x in qroots])
        print('  H roots :', [float(x) for x in Hroots])
        print('  all G>0 :', all(g > 0 for g in Gvals))


def far_lag_signs(j=2, n=8):
    """Reproduce the far-lag residue sign failures (termwise positivity false)."""
    O, H = rows(j, n)
    q = O[n - j]
    qroots = [sp.re(r) for r in sp.nroots(q, n=40) if abs(sp.im(r)) < 1e-25]
    print(f'=== far-lag residue signs, j={j} n={n} q=O_{n-j} ===')
    for k in range(0, n - j):
        lag = n - j - k
        vals = []
        for mu in qroots:
            qp = sp.diff(q, t).subs(t, mu)
            vals.append(sp.N(O[k].subs(t, mu) / qp, 12))
        has_pos = any(v > 0 for v in vals)
        has_neg = any(v < 0 for v in vals)
        print(f'  k={k} lag={lag}: '
              f'{"POS" if has_pos else ""}{" NEG" if has_neg else ""} '
              f'{[float(v) for v in vals]}')


if __name__ == '__main__':
    check_identities()
    report()
    far_lag_signs()
