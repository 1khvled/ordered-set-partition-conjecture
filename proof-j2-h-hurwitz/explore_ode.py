"""Search for a second-order ODE for P_n (evidence/exact fit)."""
import sympy as sp
from common import Family, t
F = Family(24)

# ansatz: A(t) P'' + B(t) P' + C(t) P = 0, degrees <= 2 in t,
# coefficients rational functions of n fitted at several n.
N = sp.symbols('N')
degA=degB=degC=2
unks=[]
for name,d in [('A',degA),('B',degB),('C',degC)]:
    for k in range(d+1):
        unks.append(sp.Symbol(f'{name}{k}'))

# Fit coefficients as polynomials in n of degree <=2 for each t-power
def fit_ode(parity):
    ns = [n for n in range(6,17) if n%2==parity][:6]
    eqs=[]
    for ni in ns:
        P=F.Qpoly(ni)
        expr = sum(sp.Symbol(f'A{k}')*t**k for k in range(degA+1))*P.diff().diff().as_expr()              + sum(sp.Symbol(f'B{k}')*t**k for k in range(degB+1))*P.diff().as_expr()              + sum(sp.Symbol(f'C{k}')*t**k for k in range(degC+1))*P.as_expr()
        pp=sp.Poly(sp.expand(expr),t)
        eqs.extend([c for c in pp.all_coeffs()])
    sol=sp.solve(eqs, unks, dict=True)
    return sol
for par in (0,1):
    sol=fit_ode(par)
    print('parity',par,'solutions',len(sol))
    if sol:
        for k,v in sorted(sol[0].items(), key=lambda kv:str(kv[0])):
            print(' ',k,':',sp.factor(v))
        A=sum(sol[0][sp.Symbol(f'A{k}')]*t**k for k in range(degA+1))
        B=sum(sol[0][sp.Symbol(f'B{k}')]*t**k for k in range(degB+1))
        C=sum(sol[0][sp.Symbol(f'C{k}')]*t**k for k in range(degC+1))
        print(' ODE A=',sp.factor(A),'B=',sp.factor(B),'C=',sp.factor(C))
