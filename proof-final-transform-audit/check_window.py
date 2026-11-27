from sympy import Rational, factorial, binomial, expand, Poly, symbols, nroots, re, im, Abs
from sympy import sign, N as Neval
import itertools, random

t=symbols('t')

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

def real_roots(p):
    p=expand(p)
    if p==0: return None
    Pp=Poly(p,t)
    deg=Pp.degree()
    if deg==0: return []
    # use nroots with high precision and verify
    roots=Pp.nroots(n=80)
    rr=[]
    for r in roots:
        if abs(im(r)) > 1e-30:
            return None
        rr.append(float(re(r)))
    rr.sort()
    # check distinctness approximately
    return rr

def interlaces(p,q):
    # standard: roots of p,q alternate, |deg p - deg q|<=1
    rp=real_roots(p); rq=real_roots(q)
    if rp is None or rq is None: return False, 'nonreal'
    dp=Poly(expand(p),t).degree() if expand(p)!=0 else -1
    dq=Poly(expand(q),t).degree() if expand(q)!=0 else -1
    if abs(dp-dq)>1: return False, f'degdiff {dp},{dq}'
    # merge
    # make p the larger or equal degree
    if dp<dq: dp,dq,rp,rq = dq,dp,rq,rp
    if dp==0: return True, 'const'
    # check alternation: between consecutive roots of p there is exactly one of q
    # for equal degree, also check
    ok=True
    if dp==dq:
        # interlace: r_p1 <= r_q1 <= r_p2 ...
        # check both orders
        def alt(a,b):
            # a has same length as b
            for i in range(len(a)):
                if b[i] < a[i]-1e-20: return False
                if i+1<len(a) and b[i] > a[i+1]+1e-20: return False
            return True
        ok = alt(rp,rq) or alt(rq,rp)
    else: # dp==dq+1
        for i in range(len(rp)-1):
            # need one rq in [rp[i],rp[i+1]]
            cnt=sum(1 for x in rq if x>=rp[i]-1e-20 and x<=rp[i+1]+1e-20)
            if cnt!=1: ok=False
    return ok, ('ok' if ok else f'no {rp} {rq}')

random.seed(0)
for j in [2,3,4]:
    N=14
    Hrows,O=compute_H_O(j,N)
    print(f'j={j}')
    for n in range(j, min(N+1,13)):
        # window F_n = {H_n} union {C(n,r)*O_{n-r} : 1<=r<=j-1}
        fam=[('H',n,expand(Hrows[n]))]
        for r in range(1,j):
            fam.append((f'O{n-r}',n-r,expand(binomial(n,r)*O[n-r])))
        # check real-rootedness
        for name,idx,p in fam:
            rr=real_roots(p)
            if rr is None:
                print(f'  n={n} {name} NOT REAL ROOTED')
        # pairwise interlacing
        bad=[]
        for a,b in itertools.combinations(fam,2):
            ok,msg=interlaces(a[2],b[2])
            if not ok:
                bad.append((a[0],b[0],msg))
        if bad:
            print(f'  n={n} pairwise failures: {bad}')
        # random positive combinations
        for _ in range(60):
            coeffs=[Rational(random.randint(1,5),random.randint(1,3)) for _ in fam]
            comb=expand(sum(c*p for c,(_,_,p) in zip(coeffs,fam)))
            rr=real_roots(comb)
            if rr is None:
                print(f'  n={n} COMB NOT REAL ROOTED')
                break
    print('ok')
