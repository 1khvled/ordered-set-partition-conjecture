"""Look for identities for L[R] = x^2 R''+(a+2)xR'+a^2 R and relation to lower polys."""
import sympy as sp
from ops import a, x, R

av = sp.Symbol('av', positive=True)
for n in range(4, 10):
    Rn = R(n, av)
    L = sp.expand(x**2*sp.diff(Rn,x,2) + (av+2)*x*sp.diff(Rn,x) + av**2*Rn)
    print(f"n={n}: L[R]/(leading) =", sp.factor(L / sp.LC(sp.Poly(Rn,x))))
    # try divisibility by lower Qs
    for k in range(2, n):
        for b in [av, av+1, av+2]:
            Qb = R(k, b)
            q, r = sp.div(sp.expand(L), sp.expand(Qb), x)
            if r == 0:
                print(f"   divisible by Q({k},{b})")
    print()
