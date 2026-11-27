"""Verify identity B2 = x^2[AC(c+xc') - W(2+xc)] and inspect signs."""
import sympy as sp
from explore import T
xs = sp.symbols('x')

for r in range(5, 11):
    A = T[r-1].as_expr(); C = T[r-2].as_expr()
    Ap, Cp = sp.diff(A, xs), sp.diff(C, xs)
    B2 = xs**2*(3*A*Cp - 2*Ap*C) + xs**3*(A*sp.diff(Cp, xs) - Cp*Ap)
    W = Ap*C - A*Cp
    c = Cp/C
    cp = sp.diff(c, xs)
    rhs = xs**2*(A*C*(c+xs*cp) - W*(2+xs*c))
    assert sp.simplify(sp.together(B2-rhs)) == 0, r
    pts = [-100, -3, -1, -0.1, -0.01, 0, 0.1, 1]
    print('r=%d' % r)
    print('  sign AC:', [1 if float(A.subs(xs,v)*C.subs(xs,v))>0 else -1 for v in pts])
    print('  sign W :', [1 if float(W.subs(xs,v))>0 else -1 for v in pts])
    print('  c+xc\':', ['%.3g' % float((c+xs*cp).subs(xs,v)) for v in pts])
    print('  2+xc  :', ['%.3g' % float((2+xs*c).subs(xs,v)) for v in pts])
    print('  B2    :', ['%.3g' % float(B2.subs(xs,v)) for v in pts])
