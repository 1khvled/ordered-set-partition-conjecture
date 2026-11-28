import sympy as sp
t,a=sp.symbols('t a')
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

def E(f): return sp.expand(f+t*sp.diff(f,t))
def T_op(f): return sp.expand(sp.diff(t*sp.diff(f,t),t))
def W(f,g): return sp.expand(sp.diff(f,t)*g - f*sp.diff(g,t))

# Try to express W(Pn,Tn) in terms of S=P_{n-1}, K=P_{n-2,a+1}, TS=T(S), TK terms
for n in [6,7,8]:
    c=(n-1)*a
    S=Q(n-1,a); K=Q(n-2,a+1)
    Pn=sp.expand(E(S)+c*t*K)
    Tn=sp.expand((E+1)*T_op(S) if False else T_op(Pn))
    # use actual
    Pn2=Q(n,a); Tn2=T_op(Pn2)
    print(n, sp.expand(Pn-Pn2)==0, sp.expand(Tn2 - ((E(T_op(S))+T_op(S)) + c*E(E(K))))==0)
    H=W(Pn2,Tn2)
    # try decomposition: H = A1*W(S,TS) + ...
    TS=T_op(S)
    # unknown coeffs: try to find linear combination
    # just print H at a=1
    H1=sp.expand(H.subs(a,1))
    print(' H1 deg', sp.degree(H1,t))
