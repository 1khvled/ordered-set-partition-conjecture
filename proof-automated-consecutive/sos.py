import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs
import sympy as sp

t = sp.Symbol('t')

def sos_decomp(Wexpr):
    Wp=sp.Poly(Wexpr,t)
    D=Wp.degree()
    m=D//2
    # A deg m, B deg m-1
    a=sp.symbols('a0:%d'%(m+1))
    b=sp.symbols('b0:%d'%m)
    A=sum(a[i]*t**i for i in range(m+1))
    B=sum(b[i]*t**i for i in range(m))
    expr=sp.expand(A*A + B*B - Wexpr)
    eqs=sp.Poly(expr,t).all_coeffs()
    sol=sp.solve(eqs, list(a)+list(b), dict=True)
    return sol, A, B

O={}
for n,c in ordered_rows(2,7):
    O[n]=poly_from_coeffs(c)
n=5
W=sp.expand(O[n].as_expr()*sp.diff(O[n+1].as_expr(),t)-sp.diff(O[n].as_expr(),t)*O[n+1].as_expr())
sol,A,B=sos_decomp(W)
print("sol", sol)
if sol:
    s=sol[0]
    A2=A.subs(s); B2=B.subs(s)
    print("A=",sp.expand(A2))
    print("B=",sp.expand(B2))