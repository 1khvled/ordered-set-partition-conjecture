from sympy import *
t=symbols('t')

def compute_H_O(j, N):
    def series_mul(a,b,N):
        res=[0]*(N+1)
        for i,ai in enumerate(a):
            if ai==0: continue
            for k,bk in enumerate(b):
                if bk==0: continue
                ik=i+k
                if ik>N: break
                res[ik]=res[ik]+ai*bk
        return res
    def series_inv(a,N):
        inv0 = Rational(1)/a[0]
        res=[0]*(N+1)
        res[0]=inv0
        for n in range(1,N+1):
            s=sum(a[k]*res[n-k] for k in range(1,n+1) if k<len(a) and a[k]!=0 and res[n-k]!=0)
            res[n]= -inv0*s
        return res
    A=[Rational(0)]*(N+1)
    for r in range(1,j):
        if r<=N: A[r]=Rational(1,factorial(r))
    one_minus_A=[Rational(1)] + [(-Rational(1,factorial(r)) if 1<=r<j else Rational(0)) for r in range(1,N+1)]
    D=series_inv(one_minus_A,N)
    C=[Rational(1,factorial(r)) if r>=j else Rational(0) for r in range(N+1)]
    g=series_mul(C,D,N)
    H=[Rational(0)]*(N+1); H[0]=Rational(1)
    for n in range(1,N+1):
        s=Rational(0)
        for k in range(1,n+1):
            if g[k]!=0 and H[n-k]!=0: s=s+g[k]*H[n-k]
        H[n]=expand(t*s)
    DH=series_mul(D,H,N)
    O=[expand(factorial(n)*DH[n]) for n in range(N+1)]
    Hrows=[expand(factorial(n)*H[n]) for n in range(N+1)]
    return Hrows,O

def rr(p):
    p=expand(p)
    if p==0: return None
    Pp=Poly(p,t)
    if Pp.degree()<=0: return []
    rs=Pp.nroots(n=80)
    out=[]
    for r in rs:
        if abs(im(r))>1e-25: return None
        out.append(float(re(r)))
    out.sort()
    return out

def W(p,q):
    p=expand(p); q=expand(q)
    return expand(p*diff(q,t)-diff(p,t)*q)

for j in [2,3,4]:
    Hrows,O=compute_H_O(j,14)
    print(f'\n=== j={j} ===')
    for n in range(j,11):
        Hn=expand(Hrows[n])
        print(f'n={n} degH={Poly(Hn,t).degree()} Hroots={rr(Hn)}')
        for r in range(1,j):
            Or=expand(binomial(n,r)*O[n-r])
            w=W(Hn,Or)
            # check sign of w at some points? or overall sign
            # w is a polynomial; if constant sign, interlacing possible
            wr=rr(w)
            # also check interlacing order
            rH=rr(Hn); rO=rr(Or)
            print(f'  vs B{r} deg={Poly(Or,t).degree()} roots={rO}')
            print(f'    W deg={Poly(w,t).degree() if w!=0 else -1} W_real_roots={wr}')
            # check sign of W at sample points between -10 and 0
            if w!=0:
                pts=[-10,-5,-2,-1,-0.5,-0.1,-0.01,0.0]
                vals=[float(N(w.subs(t,p))) for p in pts]
                signs=[( '+' if v>1e-12 else '-' if v<-1e-12 else '0') for v in vals]
                print(f'    Wsigns={signs}')
