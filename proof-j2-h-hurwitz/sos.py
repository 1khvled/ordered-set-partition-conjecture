"""Try SOS ansatz H = A^2 + t B^2 for actual H_n."""
import sympy as sp
from common import Family, t
F=Family(16)
for n in (6,8,10):
    H=F.H(n).as_expr(); d=sp.degree(H,t)
    # try deg A = d//2, deg B = d//2 -1 when d even (so t B^2 degree d-1)
    da=d//2; db=d//2-1
    aco=sp.symbols('a0:'+str(da+1)); bco=sp.symbols('b0:'+str(db+1))
    A=sum(aco[i]*t**i for i in range(da+1)); B=sum(bco[i]*t**i for i in range(db+1))
    eq=sp.Poly(sp.expand(H-A**2-t*B**2),t).all_coeffs()
    sol=sp.solve(eq, list(aco)+list(bco), dict=True)
    print('n',n,'deg',d,'solutions',len(sol))
    if sol:
        print('  ',sol[:2])
