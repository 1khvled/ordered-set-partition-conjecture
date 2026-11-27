from sympy import Rational, factorial, binomial, expand, Poly, symbols, nroots, re, im, gcd, LC
from itertools import combinations
import random, json, sys

t = symbols('t')
random.seed(1)

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

def real_roots(p):
    p=expand(p)
    if p==0: return None
    Pp=Poly(p,t)
    deg=Pp.degree()
    if deg<=0: return []
    roots=Pp.nroots(n=80)
    rr=[]
    for r in roots:
        if abs(im(r)) > 1e-25:
            return None
        rr.append(float(re(r)))
    rr.sort()
    # dedup very close
    return rr

def poly_deg(p):
    p=expand(p)
    if p==0: return -1
    return Poly(p,t).degree()

def leading_sign(p):
    p=expand(p)
    if p==0: return 0
    Pp=Poly(p,t)
    lc=Pp.LC()
    return 1 if lc>0 else -1 if lc<0 else 0

def interlaces(p,q):
    # returns (ok, msg, details)
    rp=real_roots(p); rq=real_roots(q)
    if rp is None or rq is None: return False,'nonreal',{}
    dp=poly_deg(p); dq=poly_deg(q)
    if abs(dp-dq)>1: return False,f'deg_jump {dp},{dq}',{}
    # normalize to leading coeff positive for sign checks
    def norm(rts,deg): return rts
    if dp<dq or (dp==dq and leading_sign(p)<0):
        # swap so that p has >= deg and prefer positive lc on larger
        pass
    # standard interlacing: smaller/eq degree interlaces larger
    if dp==dq:
        # same degree: p << q or q << p or equivalent? use common interlacer definition via Wronskian sign? 
        # For real-rooted same degree, interlace means roots alternate: p1<=q1<=p2<=...
        def alt(a,b):
            # a,b same length
            for i in range(len(a)):
                if b[i] < a[i]-1e-18: return False
                if i+1<len(a) and b[i] > a[i+1]+1e-18: return False
            return True
        ok = alt(rp,rq) or alt(rq,rp)
        order = 'p<=q' if alt(rp,rq) else ('q<=p' if alt(rq,rp) else 'none')
        return ok, order, {'rp':rp,'rq':rq}
    else:
        # different degree by 1: larger is p if dp>dq
        if dp<dq: dp,dq,rp,rq,p,q = dq,dp,rq,rp,q,p
        # dq == dp-1, check rq between rp
        ok=True
        for i in range(len(rp)-1):
            cnt=sum(1 for x in rq if x>=rp[i]-1e-18 and x<=rp[i+1]+1e-18)
            if cnt!=1: ok=False; break
        # also allow boundaries
        return ok, ('ok' if ok else 'no'), {'rp':rp,'rq':rq}

def has_common_interlacer(polys, samples=200):
    # necessary/sufficient for equal-degree family: all positive combos real-rooted.
    # for mixed degree, check all positive combos real-rooted as proxy (necessary for common interlacer of maxdeg+1?)
    # We'll check random positive combinations and also pairwise.
    for _ in range(samples):
        coeffs=[Rational(random.randint(1,7),random.randint(1,4)) for _ in polys]
        comb=expand(sum(c*p for c,p in zip(coeffs,polys)))
        if real_roots(comb) is None:
            return False
    return True

def wron_sign(p,q):
    # W[f,g] = f g' - f' g ; for interlacing, constant sign on R is related
    p=expand(p); q=expand(q)
    W=expand(p*Poly(q,t).diff() - Poly(p,t).diff()*q)
    return W

for j in [2,3,4,5]:
    N=16
    Hrows,O=compute_H_O(j,N)
    print(f'\n=== j={j} ===')
    for n in range(j, min(N,14)):
        fam=[]
        fam.append(('H',n,expand(Hrows[n])))
        for r in range(1,j):
            idx=n-r
            fam.append((f'B{r}=C({n},{r})O{idx}', idx, expand(binomial(n,r)*O[idx])))
        # degrees
        degs=[(name,poly_deg(p)) for name,_,p in fam]
        # pairwise
        bad=[]
        for (n1,i1,p1),(n2,i2,p2) in combinations(fam,2):
            ok,msg,_=interlaces(p1,p2)
            if not ok:
                bad.append((n1,n2,msg))
        ci = has_common_interlacer([p for _,_,p in fam], samples=120)
        # check O_n alignment into W_{n+1}
        On=expand(O[n])
        On_deg=poly_deg(On)
        # next window tail: C(n+1,1)O_n, ..., C(n+1,j-1)O_{n+1-j+?} wait
        # W_{n+1} = {H_{n+1}} U {C(n+1,r) O_{n+1-r} : r=1..j-1}
        # = {H_{n+1}, C(n+1,1)O_n, C(n+1,2)O_{n-1}, ..., C(n+1,j-1)O_{n+2-j}}
        next_tail=[]
        for r in range(1,j):
            idx=n+1-r
            next_tail.append((f'C({n+1},{r})O{idx}', idx, expand(binomial(n+1,r)*O[idx])))
        Hnext=('H_{n+1}', n+1, expand(Hrows[n+1]))
        # check On interlaces carried-over tail elements O_{n-1}..O_{n+2-j}
        carry_bad=[]
        for name,idx,p in next_tail:
            if idx==n: continue # that's O_n itself
            ok,msg,_=interlaces(On,p)
            if not ok:
                carry_bad.append((name,msg))
        # check Hnext interlaces On and carried tail
        hbad=[]
        for name,idx,p in [( 'O_n',n,On)] + [x for x in next_tail if x[1]!=n]:
            ok,msg,_=interlaces(Hnext[2],p)
            if not ok:
                hbad.append((name,msg))
        print(f'n={n}: degs={degs} pair_bad={bad} CI_samples={ci} carry_bad={carry_bad} Hnext_bad={hbad}')
