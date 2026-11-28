import sympy as sp
x,t,a = sp.symbols('x t a')
N=16
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

for n in range(6,13):
    for av in [1,2]:
        P=Q(n,av); TT=T_op(P); H=W(P,TT)
        signs = [sp.sign(c) for c in sp.Poly(H,t).all_coeffs()]
        print(n,av,'t-signs',signs)
        # check H(t)>0 for t<0 by counting roots in (-oo,0)
        nr = sp.Poly(H,t).count_roots(-sp.oo,0)
        if nr: print('  BAD',nr)
