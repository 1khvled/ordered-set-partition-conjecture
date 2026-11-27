"""Exact verification of W_r = B1 + (r-1) B2 and -B2 = -(r-2)x^2 W(D_C, x D_E)."""
import sympy as sp
from explore import T
xs = sp.symbols('x')
for r in range(5, 13):
    A = T[r-1].as_expr(); C = T[r-2].as_expr()
    Ap, Cp = sp.diff(A,xs), sp.diff(C,xs)
    B1 = xs**2*(A*sp.diff(Ap,xs) - Ap**2) - A**2
    B2 = xs**2*(3*A*Cp - 2*Ap*C) + xs**3*(A*sp.diff(Cp,xs) - Cp*Ap)
    F = T[r].as_expr(); G = xs*A
    W = sp.expand(sp.diff(F,xs)*G - F*sp.diff(G,xs))
    assert sp.expand(W - (B1 + (r-1)*B2)) == 0, ('decomp', r)
    Dc = 2*C + xs*Cp
    E = T[r-3].as_expr()
    De = 2*E + xs*sp.diff(E,xs)
    Wde = sp.expand(sp.diff(Dc,xs)*(xs*De) - Dc*sp.diff(xs*De,xs))
    Wcd = sp.expand(sp.diff(C,xs)*Dc - C*sp.diff(Dc,xs))
    # find alpha, beta with -B2 = xs^2 (alpha Wcd + beta Wde)
    alpha, beta = sp.symbols('alpha beta')
    sol = sp.solve(sp.Poly(sp.expand(-B2 - xs**2*(alpha*Wcd + beta*Wde)), xs).coeffs(), [alpha,beta], dict=True)
    print('r=%d decomp W_r=B1+(r-1)B2 verified; -B2 = x^2( a W(C,Dc) + b W(Dc,xDe)):' % r, sol)
    # sign of -(C/Dc)'
    R2 = -sp.diff(C/Dc, xs)
    vals = [float(R2.subs(xs,v)) for v in [-100,-3,-1,-0.1,-0.01,0.01,1]]
    print('   -(C/Dc)\' signs:', [1 if v>0 else -1 for v in vals])
