import sympy as sp
t,a=sp.symbols('t a')
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
    c=(n-1)*av
    S=Q(n-1,av); K=Q(n-2,av+1)
    TS=T_op(S)
    Pn=sp.expand(E(S)+c*t*K)
    Tn=sp.expand((E(TS)+TS)+c*E(E(K)))
    H=W(Pn,Tn)
    # expand as polynomial in c symbolically by rebuilding
    cS=sp.Symbol('cS')
    PnS=sp.expand(E(S)+cS*t*K)
    TnS=sp.expand((E(TS)+TS)+cS*E(E(K)))
    HS=W(PnS,TnS)
    HS=sp.expand(HS)
    h0=sp.expand(HS.subs(cS,0))
    # h2 coeff of cS^2
    h2=sp.expand(sp.diff(HS,cS,2)/2)
    h1=sp.expand(sp.diff(HS,cS)-2*cS*h2).subs(cS,0)
    # check signs of coeffs for each part at av=1
    def allpos(p):
        cs=sp.Poly(p,t).all_coeffs()
        return all(c>=0 for c in cs), min(cs)
    p0,m0=allpos(h0)
    p1,m1=allpos(h1)
    p2,m2=allpos(h2)
    print(n,'h0',p0,m0,'h1',p1,m1,'h2',p2,m2)
