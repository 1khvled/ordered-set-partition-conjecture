"""Exploration of threshold-j relations for O_{n,j}. Exact integer arithmetic."""
import sys
from math import comb
import sympy as sp


def O(nmax, j):
    t = None
    x = sp.symbols('t')
    O = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        val = 0
        for r in range(1, min(j - 1, n) + 1):
            val += comb(n, r) * O[n - r]
        for r in range(j, n + 1):
            val += comb(n, r) * x * O[n - r]
        O.append(sp.expand(val))
    return O, x


for j in range(1, 6):
    arr, t = O(9, j)
    print('j=', j)
    for n, p in enumerate(arr):
        print(' ', n, sp.factor(p) if p != 0 else 0)
