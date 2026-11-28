import sympy as sp
x,t,a = sp.symbols('x t a')

N=20
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

def E(f): return sp.expand(f + t*sp.diff(f,t))
def T_op(f): return sp.expand(sp.diff(t*sp.diff(f,t),t))
def W(f,g): return sp.expand(sp.diff(f,t)*g - f*sp.diff(g,t))

# Check target: W_t(P, T_P) > 0 for t<0? i.e. W_x(R,C)>0 for x>0
for n in range(4,13):
    for av in [1,2,3]:
        P = Q(n,av)
        TT = T_op(P)
        H = W(P,TT)  # W_t(P,T) 
        Hx = sp.expand(H.subs(t,-x))
        poly = sp.Poly(Hx,x)
        nr = poly.count_roots(0, sp.oo)
        if nr!=0:
            print('BAD',n,av,nr)
print('ok small check')

for n in range(6,11):
    P = Q(n,1)
    TT = T_op(P)
    H = W(P,TT)
    Hx = sp.expand(H.subs(t,-x))
    print(n, 'deg', sp.degree(Hx,x), 'lc', sp.LC(sp.Poly(Hx,x)))
    print('  coeffs signs', [sp.sign(c) for c in sp.Poly(Hx,x).all_coeffs()])
