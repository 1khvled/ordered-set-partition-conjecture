import sympy as sp
x,t,a=sp.symbols('x t a')
N=30
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

bad=[]
for n in range(4,25):
    for av in [1,2,3,5]:
        P=Q(n,av); TT=T_op(P); H=W(P,TT)
        coeffs = sp.Poly(H,t).all_coeffs()
        if any(c<=0 for c in coeffs):
            bad.append((n,av,[str(c) for c in coeffs if c<=0][:3]))
print('bad',bad)
print('checked up to n=24,a in {1,2,3,5}')
