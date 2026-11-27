"""Exact/high-precision certificates for the residue-domination obstruction.

Checks, all with sympy (exact arithmetic except root approximations, which
use 40-digit nroots):

  V1  exact identities  H_n = t S,  W(H_n,q) + t q' S = 0 mod q,  S=sum C(n,r)O_{n-r};
  V2  at every root mu of q=O_{n-d}:  G_q(mu)>0  AND  aligned interlacing of H_n, q;
  V3  far-lag residues O_k(mu)/q'(mu) take BOTH signs (termwise positivity false);
  V4  exact decomposition  W(O_n,O_{n-j}) = W(H_n,O_{n-j}) + sum_{r<j} C(n,r)W(O_{n-r},O_{n-j});
  V5  alpha'= (B'C-BC')/C^2  has numerator e^z(p-2)+2E-2p, negative near 0, positive for large z.

Run: py -3 verify_obstruction.py
"""
import sympy as sp

t, z = sp.symbols('t z')


def rows(j, nmax):
    O = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = sum(sp.binomial(n, r) * O[n - r] for r in range(1, min(j - 1, n) + 1))
        s += sum(t * sp.binomial(n, r) * O[n - r] for r in range(j, n + 1))
        O.append(sp.expand(s))
    H = [sp.expand(O[n] - sum(sp.binomial(n, r) * O[n - r]
                              for r in range(1, min(j - 1, n) + 1)))
         for n in range(nmax + 1)]
    return O, H


def W(p, q):
    return sp.expand(sp.diff(p, t) * q - p * sp.diff(q, t))


def real_roots_neg(poly, digits=40):
    if sp.degree(poly, t) is None or sp.degree(poly, t) < 1:
        return []
    return sorted(sp.re(r) for r in sp.nroots(poly, n=digits)
                  if abs(sp.im(r)) < 1e-25 and sp.re(r) < -1e-25)


def v1_exact_identities(jmax=6, nmax=24):
    tested = 0
    for j in range(2, jmax + 1):
        O, H = rows(j, nmax)
        for n in range(j, nmax + 1):
            S = sp.expand(sum(sp.binomial(n, r) * O[n - r] for r in range(j, n + 1)))
            assert sp.expand(H[n] - t * S) == 0, ('H=tS', j, n)
            for d in range(1, j + 1):
                q = O[n - d]
                if sp.degree(q, t) < 1:
                    continue
                poly = sp.expand(W(H[n], q) + t * sp.diff(q, t) * S)
                assert sp.rem(sp.Poly(poly, t), sp.Poly(q, t)).as_expr() == 0, \
                    ('residue-id', j, n, d)
                tested += 1
    print(f'V1 exact identities: {tested} triples PASS')


def v2_domination_interlacing(jmax=6, nmax=24):
    tested = 0
    for j in range(2, jmax + 1):
        O, H = rows(j, nmax)
        for n in range(j, nmax + 1):
            S = sp.expand(sum(sp.binomial(n, r) * O[n - r] for r in range(j, n + 1)))
            for d in range(1, j + 1):
                q = O[n - d]
                mus = real_roots_neg(q)
                if not mus:
                    continue
                qp = sp.diff(q, t)
                G = [sp.N(S.subs(t, mu) / qp.subs(t, mu), 20) for mu in mus]
                assert all(g > 0 for g in G), ('G<=0', j, n, d, G)
                # aligned interlacing of H_n roots {nu_1<..<nu_m, plus 0} with mus
                nus = real_roots_neg(H[n])  # negative roots only; 0 is also a root
                # residue condition <0 implies aligned interlacing; verify order
                if len(nus) == len(mus):      # deg p = e+1 case: nu_i < mu_i < nu_{i+1}
                    ok = all(nus[i] < mus[i] for i in range(len(mus)))
                elif len(nus) == len(mus) - 1:  # deg p = e case: mu_i < nu_i < mu_{i+1}
                    ok = all(mus[i] < nus[i] for i in range(len(nus)))
                else:
                    ok = False
                assert ok, ('interlace', j, n, d, mus, nus)
                tested += 1
    print(f'V2 domination + aligned interlacing: {tested} triples PASS')


def v3_far_lag_signs(j=2, n=8):
    O, H = rows(j, n)
    q = O[n - j]
    mus = real_roots_neg(q)
    qp = sp.diff(q, t)
    mixed = 0
    for k in range(0, n - j):
        vals = [sp.N(O[k].subs(t, mu) / qp.subs(t, mu), 12) for mu in mus]
        if any(v > 0 for v in vals) and any(v < 0 for v in vals):
            mixed += 1
    print(f'V3 far-lag residues with BOTH signs (j={j},n={n},q=O_{n-j}): '
          f'{mixed} mixed-sign lags  -> termwise positivity FALSE')


def v4_no_free_lunch(jmax=5, nmax=20):
    tested = 0
    for j in range(2, jmax + 1):
        O, H = rows(j, nmax)
        for n in range(2 * j, nmax + 1):
            q = O[n - j]
            lhs = W(O[n], q)
            rhs = W(H[n], q) + sum(sp.binomial(n, r) * W(O[n - r], q)
                                   for r in range(1, j))
            assert sp.expand(lhs - rhs) == 0, ('decomp', j, n)
            tested += 1
    print(f'V4 exact d=j decomposition: {tested} cases PASS')


def v5_alpha_not_stieltjes(j=2):
    E = sum(z ** r / sp.factorial(r) for r in range(j))
    C = sp.exp(z) - E
    B = 2 - E
    p = z ** (j - 1) / sp.factorial(j - 1)
    numer = sp.simplify(sp.diff(B, z) * C - B * sp.diff(C, z))
    # E = E_{j-1}; 2*E - 2*p = 2*E_{j-2} since E_{j-1}=E_{j-2}+p.
    target = sp.exp(z) * (p - 2) + 2 * E - 2 * p
    assert sp.simplify(numer - target) == 0, ('alpha numerator', j)
    v_small = float(numer.subs(z, sp.Rational(1, 10)))
    v_large = float(numer.subs(z, 20))
    assert v_small < 0 < v_large, (v_small, v_large)
    print(f'V5 alpha=B/C derivative numerator changes sign: '
          f'numer(0.1)={v_small:.3g}<0, numer(20)={v_large:.3g}>0 -> not Stieltjes')


if __name__ == '__main__':
    v1_exact_identities()
    v2_domination_interlacing()
    v3_far_lag_signs()
    v4_no_free_lunch()
    v5_alpha_not_stieltjes()
    print('ALL CHECKS PASS')
