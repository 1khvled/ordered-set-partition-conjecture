from sympy import Rational,factorial,binomial,expand,Poly,symbols,nroots,re,im
t=symbols('t')
# reuse compute
def compute_H_O(j, N):
    def series_mul(a,b,N):
        res=[0]*(N+1)
        for i,ai in enumerate(a):
            if ai==0: continue
            for k,bk in enumerate(b):
                if bk==0: continue
                if i+k>N: break
                res[i+k]=res[i+k]+ai*bk
        return res
    def series_inv(a,N):
        inv0 = 1/a[0]
        res=[0]*(N+1)
        res[0]=inv0
        for n in range(1,N+1):
            s=sum(a[k]*res[n-k] for k in range(1,n+1) if k<len(a) and a[k]!=0 and res[n-k]!=0)
            res[n]= -inv0*s
        return res
    A=[Rational(0)]*(N+1)
    for r in range(1,j):
        if r<=N: A[r]=Rational(1,factorial(r))
    one_minus_A=[Rational(1)] + [-Rational(1,factorial(r)) if 1<=r<j else Rational(0) for r in range(1,N+1)]
    D=series_inv(one_minus_A,N)
    C=[Rational(1,factorial(r)) if r>=j else Rational(0) for r in range(N+1)]
    g=series_mul(C,D,N)
    H=[0]*(N+1); H[0]=1
    for n in range(1,N+1):
        s=0
        for k in range(1,n+1):
            if g[k]!=0 and H[n-k]!=0: s=s+g[k]*H[n-k]
        H[n]=expand(t*s)
    DH=series_mul(D,H,N)
    O=[expand(factorial(n)*DH[n]) for n in range(N+1)]
    Hrows=[expand(factorial(n)*H[n]) for n in range(N+1)]
    return Hrows,O

def rroots(p):
    Pp=Poly(expand(p),t)
    if Pp.degree()<=0: return []
    rs=[complex(r) for r in Pp.nroots(n=60)]
    rs=[float(re(r)) for r in rs if abs(im(r))<1e-25]
    rs.sort()
    return rs

for j in [2,3]:
    Hrows,O=compute_H_O(j,12)
    for n in range(j,10):
        print(f'j={j} n={n}')
        print(' H', rroots(Hrows[n]), 'deg', Poly(Hrows[n],t).degree())
        for r in range(1,j):
            p=expand(binomial(n,r)*O[n-r])
            print(f' O_{n-r}', rroots(p), 'deg', Poly(p,t).degree(), 'val0', p.subs(t,0))
