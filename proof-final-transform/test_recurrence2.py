import transform_lib as tl
import sympy as sp
t=sp.symbols('t')
for j in [2,3,4,5]:
    N=8
    O,H,D=tl.compute_O(j,N)
    ok=True
    for n in range(N+1):
        lhs=O[n].as_expr()
        rhs=sp.factorial(n)*H[n].as_expr()
        for r in range(1,j):
            if n-r>=0:
                rhs=rhs+sp.binomial(n,r)*O[n-r].as_expr()
        if sp.expand(lhs-rhs)!=0:
            ok=False; print('MISMATCH j',j,'n',n)
    print('j',j,'recurrence verified:',ok)
