import transform_lib as tl
import sympy as sp
import random
t=sp.symbols('t')
for j in [4,5]:
    N=8
    O,H,D=tl.compute_O(j,N)
    def family(n):
        F=[H[n]]
        for r in range(1,j):
            if n-r>=0:
                F.append(sp.Poly(sp.binomial(n,r)*O[n-r].as_expr(), t))
        return F
    for n in range(j, N+1):
        F=family(n)
        bad=0
        for trial in range(100):
            weights=[random.random()+0.1 for _ in F]
            combo=sum(w*F[i].as_expr() for i,w in enumerate(weights))
            if not tl.is_real_rooted(sp.expand(combo)):
                bad+=1
        if bad>0: print('j',j,'n',n,'bad',bad)
    print('j',j,'done')
