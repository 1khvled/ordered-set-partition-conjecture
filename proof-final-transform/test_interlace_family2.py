import transform_lib as tl
import sympy as sp
import random
t=sp.symbols('t')
for j in [2,4,5]:
    N=10
    O,H,D=tl.compute_O(j,N)
    def family(n):
        F=[H[n]]
        for r in range(1,j):
            if n-r>=0:
                F.append(sp.Poly(sp.binomial(n,r)*O[n-r].as_expr(), t))
        return F
    def real_rooted_expr(expr):
        return tl.is_real_rooted(expr)
    worst=0
    for n in range(j, N+1):
        F=family(n)
        bad=0
        for trial in range(300):
            weights=[random.random()+0.1 for _ in F]
            combo=sum(w*F[i].as_expr() for i,w in enumerate(weights))
            if not real_rooted_expr(sp.expand(combo)):
                bad+=1
        worst=max(worst,bad)
        if bad>0: print('j',j,'n',n,'bad',bad)
    print('j',j,'worst bad',worst)
