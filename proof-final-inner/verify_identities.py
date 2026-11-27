"""Exploration/verification only: exact checks of identities for the inner family.
Not a proof of real-rootedness.
"""
import math
from fractions import Fraction
from explore_H import inner_array


def check_recurrence(j, nmax):
    h = inner_array(j, nmax)
    # inner_array stores h[n][m]/n!; convert to the integer array h_{n,m}
    h = [[Fraction(math.factorial(n)) * h[n][m] for m in range(len(h[n]))]
         for n in range(len(h))]
    errs = []
    for n in range(1, nmax + 1):
        for m in range(1, n // j + 2):
            # h[n][m] = sum_{k=1}^{j-1} C(n,k) h[n-k][m]
            #           + sum_{q=j}^{n} C(n,q) h[n-q][m-1]
            lhs = h[n][m]
            rhs = Fraction(0)
            for k in range(1, min(j, n + 1)):
                rhs += math.comb(n, k) * h[n - k][m]
            for q in range(j, n + 1):
                rhs += math.comb(n, q) * h[n - q][m - 1]
            if lhs != rhs:
                errs.append((n, m, lhs, rhs))
    return errs


def check_poly_rec(j, nmax):
    """H_n = sum_{k=1}^{j-1} C(n,k) H_{n-k} + t sum_{q=j}^n C(n,q) H_{n-q}, for n>=j."""
    h = inner_array(j, nmax)
    h = [[Fraction(math.factorial(n)) * h[n][m] for m in range(len(h[n]))]
         for n in range(len(h))]
    errs = []
    for n in range(j, nmax + 1):
        for m in range(0, n // j + 2):
            lhs = h[n][m]
            rhs = Fraction(0)
            for k in range(1, j):
                rhs += math.comb(n, k) * h[n - k][m]
            for q in range(j, n + 1):
                if m >= 1:
                    rhs += math.comb(n, q) * h[n - q][m - 1]
            if lhs != rhs:
                errs.append((n, m, lhs, rhs))
    return errs


if __name__ == '__main__':
    for j in [2, 3, 4, 5]:
        e1 = check_recurrence(j, 18)
        e2 = check_poly_rec(j, 18)
        print('j=%d  recurrence errors: %d   poly-rec errors: %d' % (j, len(e1), len(e2)))
        if e1:
            print('  first:', e1[0])
        if e2:
            print('  first:', e2[0])
