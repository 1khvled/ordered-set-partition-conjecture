import sympy as sp
t=sp.symbols('t')
N=45
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

bad_coeff=[]
bad_real=[]
for n in range(4,41):
    for av in range(1,11):
        P=Q(n,av); TT=T_op(P); H=W(P,TT)
        coeffs=sp.Poly(H,t).all_coeffs()
        if any(c<=0 for c in coeffs):
            bad_coeff.append((n,av))
        rr=sp.Poly(H,t).count_roots(-sp.oo,sp.oo)
        if rr: bad_real.append((n,av,rr))
    if n%5==0: print('done',n,flush=True)
print('bad_coeff',bad_coeff)
print('bad_real',bad_real)
