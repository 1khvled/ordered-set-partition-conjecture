import sympy as sp
t=sp.symbols('t')
N=12
A=[[0]*(N//2+3) for _ in range(N+1)]
A[0][0]=1
for m in range(N):
    for r in range(len(A[0])-1):
        A[m+1][r]+=r*A[m][r]
        if r: A[m+1][r]+=m*A[m-1][r-1]

def rising(av,r):
    v=1
    for i in range(r): v*=av+i
    return v

def Q(n,av,var=t):
    return sp.expand(sum(rising(av,r)*A[n][r]*var**(r-1) for r in range(1,len(A[0])) if A[n][r]))
def T_op(f): return sp.expand(sp.diff(t*sp.diff(f,t),t))
def W(f,g): return sp.expand(sp.diff(f,t)*g - f*sp.diff(g,t))

for n in [6,8,10]:
    P=Q(n,1); TT=T_op(P); H=W(P,TT)
    print(n, H)
    print(' factor', sp.factor(H))
    # check if positive on R by sturm on (-oo,oo)
    print(' roots on R', sp.Poly(H,t).count_roots(-sp.oo,sp.oo))
