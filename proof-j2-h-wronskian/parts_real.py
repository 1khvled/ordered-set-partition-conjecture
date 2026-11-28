import sympy as sp
t=sp.symbols('t')
N=14
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
def E(f): return sp.expand(f+t*sp.diff(f,t))
def T_op(f): return sp.expand(sp.diff(t*sp.diff(f,t),t))
def W(f,g): return sp.expand(sp.diff(f,t)*g - f*sp.diff(g,t))

for n in range(6,12):
    av=1
    S=Q(n-1,av); K=Q(n-2,av+1)
    TS=T_op(S)
    cS=sp.Symbol('cS')
    PnS=sp.expand(E(S)+cS*t*K)
    TnS=sp.expand((E(TS)+TS)+cS*E(E(K)))
    HS=sp.expand(W(PnS,TnS))
    h0=sp.expand(HS.subs(cS,0))
    h2=sp.expand(sp.diff(HS,cS,2)/2)
    h1=sp.expand(sp.diff(HS,cS)-2*cS*h2).subs(cS,0)
    for name,h in [('h0',h0),('h1',h1),('h2',h2)]:
        rr=sp.Poly(h,t).count_roots(-sp.oo,sp.oo)
        if rr: print(n,name,'REAL ROOTS',rr)
print('ok')
