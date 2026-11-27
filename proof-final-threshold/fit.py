"""Try to fit O_{n,j} as a binomial convolution of O_{k,j-1} with a kernel f_s(t)."""
from math import comb
import sympy as sp

t = sp.symbols('t')

def O(nmax, j):
    arr = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        v = 0
        for r in range(1, min(j - 1, n) + 1):
            v += comb(n, r) * arr[n - r]
        for r in range(j, n + 1):
            v += comb(n, r) * t * arr[n - r]
        arr.append(sp.expand(v))
    return arr

for j in [2, 3]:
    A = O(12, j)
    B = O(12, j - 1)
    f = [sp.Integer(1)]
    print('j=', j)
    for n in range(1, 11):
        # solve f_n from k=0 term
        rhs = sum(comb(n, k) * f[n - k] * B[k] for k in range(1, n + 1))
        fn = sp.expand(A[n] - rhs)
        f.append(fn)
        print(' f_%d =' % n, sp.factor(fn))
    # consistency check: at n+1, k=1 uses f_n
    for n in range(1, 10):
        pred = sum(comb(n + 1, k) * f[n + 1 - k] * B[k] for k in range(0, n + 2))
        if sp.expand(pred - A[n + 1]) != 0:
            print(' INCONSISTENT at n=', n + 1, sp.expand(pred - A[n + 1]))
            break
    else:
        print(' consistent (by construction at k=0 only; cross-check below)')
