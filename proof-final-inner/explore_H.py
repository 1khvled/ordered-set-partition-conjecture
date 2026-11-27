"""Exploration only (not a proof): exact coefficient arrays for the inner family."""
import math
from fractions import Fraction


def inner_array(j, nmax, mmax=None):
    """h[n][m] = [t^m] H_n, H_n = n![z^n] 1/(1-t g_j), g_j = C_j/(1-A_j)."""
    if mmax is None:
        mmax = nmax // j + 1
    # ordinary coefficients of g_j up to nmax
    A = [Fraction(0)] * (nmax + 1)
    for r in range(1, j):
        A[r] = Fraction(1, math.factorial(r))
    D = [Fraction(0)] * (nmax + 1)
    D[0] = Fraction(1)
    for n in range(1, nmax + 1):
        D[n] = sum(A[r] * D[n - r] for r in range(1, min(j, n + 1)))
    C = [Fraction(0)] * (nmax + 1)
    for r in range(j, nmax + 1):
        C[r] = Fraction(1, math.factorial(r))
    g = [Fraction(0)] * (nmax + 1)
    for n in range(nmax + 1):
        g[n] = sum(C[q] * D[n - q] for q in range(j, n + 1))
    # H = 1/(1 - t g): recurrence H_n = sum_{k=1}^n g_k H_{n-k} t
    # as polynomials in t via coefficient arrays
    h = [[Fraction(0)] * (mmax + 1) for _ in range(nmax + 1)]
    h[0][0] = Fraction(1)
    for n in range(1, nmax + 1):
        for k in range(1, n + 1):
            if g[k] == 0:
                continue
            for m in range(1, mmax + 1):
                h[n][m] += g[k] * h[n - k][m - 1]
    return h


if __name__ == '__main__':
    for j in [2, 3, 4]:
        h = inner_array(j, 16)
        print('j =', j)
        for n in range(16):
            row = [str(h[n][m]) for m in range(1, n // j + 2) if h[n][m] != 0]
            print('  n=%2d:' % n, row)
