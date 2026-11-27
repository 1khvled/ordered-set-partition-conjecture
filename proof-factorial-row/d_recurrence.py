"""Does D_r=2F_r+xF_r' satisfy the same recurrence as F_r?"""
import sympy as sp
from explore import T
xs = sp.symbols('x')
D = {r: sp.Poly(sp.expand(2*T[r].as_expr() + xs*sp.diff(T[r].as_expr(), xs)), xs) for r in range(2, 20)}
for r in range(5, 15):
    lhs = D[r]
    rhs = sp.Poly(sp.expand(sp.diff(xs*D[r-1].as_expr(), xs) + (r-1)*sp.diff(xs**2*D[r-2].as_expr(), xs)), xs)
    print('r=%d same recurrence? %s' % (r, lhs == rhs))
    if lhs != rhs:
        print('   diff:', sp.factor((lhs-rhs).as_expr()))
