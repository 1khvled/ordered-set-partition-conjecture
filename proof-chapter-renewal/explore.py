"""Exploration only (not a certificate): chapter/renewal structure of H_n, O_n."""
import sympy as sp

t, z = sp.symbols('t z')


def rows(j, N):
    """Exact O_n, H_n via first-block recurrence."""
    O = [sp.Integer(1)]
    H = [sp.Integer(1)]
    for n in range(1, N + 1):
        On = sum(sp.binomial(n, r) * (1 if r < j else t) * O[n - r]
                 for r in range(1, n + 1))
        O.append(sp.expand(On))
        Hn = t * sum(sp.binomial(n, r) * O[n - r] for r in range(j, n + 1))
        H.append(sp.expand(Hn))
    return O, H


def W(p, q):
    return sp.expand(sp.diff(p, t) * q - p * sp.diff(q, t))


if __name__ == '__main__':
    for j in [2, 3]:
        O, H = rows(j, 4 * j)
        print('j =', j)
        for n in range(j, 3 * j + 1):
            print(' n=%d O_n=%s' % (n, O[n]))
            print('     H_n=%s' % H[n])
        print('  W(H_n,O_{n-j}) for n=j..3j:')
        for n in range(j, 3 * j + 1):
            w = W(H[n], O[n - j])
            print('   n=%d  W=%s  coeffs>=0? %s' % (
                n, sp.factor(w), all(c >= 0 for c in sp.Poly(w, t).all_coeffs())))
