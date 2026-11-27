import transform_lib as tl
import sympy as sp
t=sp.symbols('t')
for j in [2,3,4,5]:
    N=12
    O,H,D=tl.compute_O(j,N)
    bad=[]
    for n in range(N+1):
        expr=O[n].as_expr()
        if not tl.is_real_rooted(expr):
            bad.append(n)
    print('j',j,'O non-real n:',bad)
    # H
    H,D,g=tl.compute_H(j,N)
    badh=[]
    for r in range(N+1):
        if not tl.is_real_rooted(H[r].as_expr()): badh.append(r)
    print('j',j,'H non-real r:',badh)
