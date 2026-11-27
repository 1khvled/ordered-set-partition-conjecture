import transform_lib as tl
import sympy as sp
import random
t=sp.symbols('t')
j=3; N=8
O,H,D=tl.compute_O(j,N)
# For a given n, family: H_n, binom(n,r) O_{n-r}
def family(n):
    F=[H[n]]
    for r in range(1,j):
        if n-r>=0:
            F.append(sp.Poly(sp.binomial(n,r)*O[n-r].as_expr(), t))
    return F
def real_rooted_expr(expr):
    return tl.is_real_rooted(expr)
for n in range(j, N+1):
    F=family(n)
    # test many random positive combos
    bad=0
    for trial in range(200):
        weights=[random.random()+0.1 for _ in F]
        combo=sum(w*F[i].as_expr() for i,w in enumerate(weights))
        if not real_rooted_expr(sp.expand(combo)):
            bad+=1
            if bad<=3: print('n',n,'bad combo',weights)
    print('j',j,'n',n,'random combo bad:',bad,'/200')
