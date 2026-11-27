"""Exploration only: coefficientwise signs of far-lag Wronskians."""
import sympy as sp
from explore import rows, W

t = sp.symbols('t')

if __name__ == '__main__':
    for j in [2, 3, 4]:
        O, H = rows(j, 5 * j)
        n = 4 * j
        q = O[n - j]
        print('=== j=%d n=%d, q=O_%d' % (j, n, n - j))
        for r in range(j, n + 1):
            w = W(O[n - r], q)
            cs = sp.Poly(w, t).all_coeffs() if w != 0 else [0]
            print('  r=%2d lag=%2d  all coeffs >=0? %s  min=%s' % (
                r, r - j, all(c >= 0 for c in cs), min(cs)))
