"""Simplify W(P,Q) bracket in log-derivative variables a,b,ap=a',bp=b'."""
import sympy as sp
xs, a, b, ap, bp = sp.symbols('x a b ap bp')
B = (3*a + xs*a**2 + xs*ap)*(2+xs*b) - (2+xs*a)*(2+4*xs*b+xs**2*b**2+xs*bp)
print('B =', sp.factor(sp.expand(B)))
print()
# substitute ap=-Sp, bp=-Sq
Sp, Sq = sp.symbols('Sp Sq', nonnegative=True)
B2 = sp.expand(B.subs({ap:-Sp, bp:-Sq}))
print('B(Sp,Sq) =', sp.collect(B2, [Sp, Sq, b, a]))
