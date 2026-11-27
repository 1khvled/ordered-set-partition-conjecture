import transform_lib as tl
import sympy as sp
t=sp.symbols('t')
O,H,D=tl.compute_O(2,4)
print('compute_O O2', O[2].as_expr())
# actual
def O_rec(j,N):
    O=[None]*(N+1); O[0]=sp.Poly(1,t)
    for n in range(1,N+1):
        s=sp.Poly(0,t)
        for r in range(1,j):
            if r<=n: s=s+sp.Poly(sp.binomial(n,r)*O[n-r].as_expr(),t)
        for r in range(j,n+1):
            s=s+sp.Poly(t*sp.binomial(n,r)*O[n-r].as_expr(),t)
        O[n]=sp.Poly(sp.expand(s.as_expr()),t)
    return O
Oa=O_rec(2,4)
print('actual O2', Oa[2].as_expr())
