"""Look for exact identity for -W_r - F_{r-1}^2."""
import sympy as sp
from explore import T
xs = sp.symbols('x')

for r in range(4, 13):
    F = T[r].as_expr()
    G = xs * T[r - 1].as_expr()
    W = sp.expand(sp.diff(F, xs) * G - F * sp.diff(G, xs))
    D = sp.factor(sp.expand(-W - T[r-1].as_expr()**2))
    print('r=%d  -W_r - F_{r-1}^2 =' % r, D)
