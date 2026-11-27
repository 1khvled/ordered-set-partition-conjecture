"""Exact checks for proof-joint-interlacing."""
import sympy as sp

t = sp.Symbol('t')

def rows(j, nmax):
    O = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = sp.Integer(0)
        for r in range(1, min(j - 1, n) + 1):
            s += sp.binomial(n, r) * O[n - r]
        for r in range(j, n + 1):
            s += t * sp.binomial(n, r) * O[n - r]
        O.append(sp.expand(s))
    H = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = sp.Integer(0)
        for r in range(j, n + 1):
            s += sp.binomial(n, r) * O[n - r]
        H.append(sp.expand(t * s))
    return O, H

def W(p, q):
    return sp.expand(sp.diff(p, t) * q - p * sp.diff(q, t))

def main():
    # 1) common recurrence
    rec_fail = 0
    rec_total = 0
    for j in range(2, 7):
        nmax = 4 * j
        O, H = rows(j, nmax)
        for n in range(j, nmax + 1):
            for X, name in [(O, 'O'), (H, 'H')]:
                rhs = sp.Integer(0)
                for k in range(1, j):
                    rhs += sp.binomial(n, k) * X[n - k]
                for r in range(j, n + 1):
                    rhs += t * sp.binomial(n, r) * X[n - r]
                rhs = sp.expand(rhs)
                rec_total += 1
                if sp.expand(X[n] - rhs) != 0:
                    rec_fail += 1
                    print('REC FAIL', j, n, name)
    print('recurrence checks', rec_total, 'fail', rec_fail)

    # 2) residue identity and decomposition
    res_fail = 0
    res_total = 0
    dec_fail = 0
    dec_total = 0
    for j in range(2, 7):
        nmax = 4 * j
        O, H = rows(j, nmax)
        for n in range(j, nmax + 1):
            for d in range(1, j + 1):
                m = n - d
                if m < 0:
                    continue
                q = O[m]
                if q == 0:
                    continue
                # residue identity: W(H,q)(mu) = -H(mu)q'(mu) at roots of q
                # check polynomial identity: W(H,q) + H*q' is divisible? no, better check S_n
                Sn = sp.Integer(0)
                for k in range(0, n - j + 1):
                    Sn += sp.binomial(n, k) * O[k]
                # H_n = t*S_n
                res_total += 1
                if sp.expand(H[n] - sp.expand(t * Sn)) != 0:
                    res_fail += 1
                    print('H=tS FAIL', j, n, d)
                # decomposition of W(O_n,q_d)
                lhs = W(O[n], q)
                rhs = W(H[n], q)
                for r in range(1, d):
                    rhs += sp.binomial(n, r) * W(O[n - r], q)
                for r in range(d + 1, j):
                    rhs -= sp.binomial(n, r) * W(q, O[n - r])
                rhs = sp.expand(rhs)
                dec_total += 1
                if sp.expand(lhs - rhs) != 0:
                    dec_fail += 1
                    print('DEC FAIL', j, n, d)
                # root residue sign check: S_n(mu)/q'(mu)>0
                qpoly = sp.Poly(q, t)
                if qpoly.degree() >= 1:
                    for mu in sp.real_roots(qpoly):
                        qp = sp.diff(q, t)
                        val = sp.simplify(Sn.subs(t, mu) / qp.subs(t, mu))
                        res_total += 1
                        # val should be > 0
                        if not (sp.N(val) > 0):
                            res_fail += 1
                            print('RES SIGN FAIL', j, n, d, mu)
    print('res checks', res_total, 'fail', res_fail)
    print('decomp checks', dec_total, 'fail', dec_fail)

if __name__ == '__main__':
    main()
