"""Search identity: W(Dc, xDe) = self(E) + k*x^2*W(UE, xUH) where C=(xE)'+m(x^2H)'."""
import sympy as sp
xs, m, k = sp.symbols('x m k')
E = sp.Function('E')(xs)
H = sp.Function('H')(xs)
U = lambda f: 2*f + xs*sp.diff(f, xs)
C = sp.diff(xs*E, xs) + m*sp.diff(xs**2*H, xs)
Dc = U(C)
De = U(E)
W1 = sp.expand(sp.diff(Dc,xs)*(xs*De) - Dc*sp.diff(xs*De,xs))
W2 = sp.expand(sp.diff(U(E),xs)*(xs*U(H)) - U(E)*sp.diff(xs*U(H),xs))
# search for lambda (rational in x) such that W1 - m*lambda*W2 has no H
lam = sp.symbols('lam')
expr = sp.expand(W1 - m*lam*W2)
for dv in [3,2,1,0]:
    coef = sp.expand(expr).coeff(sp.diff(H,xs,dv) if dv else H)
    print('H^%d coeff:' % dv, sp.factor(coef))

