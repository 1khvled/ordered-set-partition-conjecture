"""Exploration only: coefficient signs of the product-rule kernel p q + t W(p,q)."""
import sympy as sp
from explore import rows, W

t = sp.symbols('t')

if __name__ == '__main__':
    for j in [2, 3]:
        O, H = rows(j, 5 * j)
        n = 4 * j
        q = O[n - j]
        print('=== j=%d  q=O_%d' % (j, n - j))
        for r in range(j, n + 1):
            p = O[n - r]
            expr = sp.expand(p * q + t * W(p, q))
            cs = sp.Poly(expr, t).all_coeffs()
            print('  r=%2d  pq+tW coeffs>=0? %s  min=%s' % (
                r, all(c >= 0 for c in cs), min(cs)))
