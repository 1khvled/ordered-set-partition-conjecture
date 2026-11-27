from sympy import *
import random, itertools
t=symbols('t')
random.seed(2)

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

def is_real_rooted(p):
    p=expand(p)
    if p==0: return True
    Pp=Poly(p,t)
    if Pp.degree()<=0: return True
    rs=Pp.nroots(n=60)
    return all(abs(im(r))<=1e-25 for r in rs)

for j in [2,3,4,5,6]:
    N=20
    Hrows,O=compute_H_O(j,N)
    for n in range(j, N):
        fam=[expand(Hrows[n])]+[expand(binomial(n,r)*O[n-r]) for r in range(1,j)]
        # pairwise Wronskian constant sign
        ok=True
        for a,b in itertools.combinations(fam,2):
            w=expand(a*diff(b,t)-diff(a,t)*b)
            if w==0: continue
            Pp=Poly(w,t)
            if Pp.degree()>0:
                rs=Pp.nroots(n=60)
                if any(abs(im(r))<=1e-25 for r in rs):
                    # real root means sign change possible
                    print('WRONSKIAN_REAL_ROOT',j,n)
                    ok=False
        # random combos
        for _ in range(80):
            coeffs=[Rational(random.randint(1,9),random.randint(1,5)) for _ in fam]
            comb=expand(sum(c*p for c,p in zip(coeffs,fam)))
            if not is_real_rooted(comb):
                print('COMBO_FAIL',j,n)
                ok=False
                break
        # propagation
        On=expand(O[n])
        Wn1=[expand(Hrows[n+1])]+[expand(binomial(n+1,r)*O[n+1-r]) for r in range(1,j)]
        for a,b in itertools.combinations(Wn1,2):
            w=expand(a*diff(b,t)-diff(a,t)*b)
            if w==0: continue
            Pp=Poly(w,t)
            if Pp.degree()>0:
                rs=Pp.nroots(n=60)
                if any(abs(im(r))<=1e-25 for r in rs):
                    print('PROP_W_REAL',j,n+1)
        if not ok:
            print('bad',j,n)
    print('done',j)
