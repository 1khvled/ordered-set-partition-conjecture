"""Exact verification of Theorem B3 (truncated inner renewal) and (B1)."""
import sympy as sp
from explore import rows, W

t = sp.symbols('t')


def d_coeffs(j, N):
    """d_r = r![z^r] 1/(1-A_j(z)) for r <= N."""
    A = sum(z ** r / sp.factorial(r) for r in range(1, j)) if False else None
    zz = sp.symbols('zz')
    A = sum(zz ** r / sp.factorial(r) for r in range(1, j))
    ser = sp.series(1 / (1 - A), zz, 0, N + 1).removeO()
    return [sp.expand(sp.diff(ser, zz, r).subs(zz, 0)) for r in range(N + 1)]


z = sp.symbols('z')

if __name__ == '__main__':
    for j in [2, 3, 4, 5]:
        O, H = rows(j, 4 * j)
        d = d_coeffs(j, 4 * j)
        for n in range(j, 4 * j + 1):
            lhs = sp.expand(H[n] / t)
            rhs = 0
            for s in range(0, n - j + 1):
                ker = sum(sp.binomial(n - s, l) * d[l]
                          for l in range(0, n - j - s + 1))
                rhs += sp.binomial(n, s) * H[s] * ker
            assert sp.expand(lhs - rhs) == 0, (j, n, sp.expand(lhs - rhs))
        # check (B1) for d=1..j
        for n in range(j, 4 * j + 1):
            for dd in range(1, j + 1):
                q = O[n - dd]
                lhs = W(H[n], q)
                rhs = sum(sp.binomial(n, r) * (O[n - r] * q + t * W(O[n - r], q))
                          for r in range(j, n + 1))
                assert sp.expand(lhs - rhs) == 0, (j, n, dd)
        print('j=%d: B3 and B1 verified for n <= %d' % (j, 4 * j))
