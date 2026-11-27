import sys, math
sys.path.insert(0, r"C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-inner")
from explore_H import inner_array
from fractions import Fraction

def check(j,nmax):
    h = inner_array(j, nmax)
    # h[n][m] as Fraction (already [t^m] H_n? let's verify by printing)
    # We'll treat h[n] as coeff list
    errs=[]
    for n in range(1, nmax):
        deg_prev = (n)//j
        # compute H_{n+1} via recurrence using T,Q
        # Represent polynomials as lists low->high of Fractions
        def get(N):
            row=h[N]
            d=N//j
            return [row[m] for m in range(d+1)]
        # easier: work coefficients
        def deriv(p):
            return [ (k+1)*p[k+1] for k in range(len(p)-1) ] if len(p)>1 else [Fraction(0)]
        def mul_t(p):
            return [Fraction(0)]+p
        def add(a,b):
            n=max(len(a),len(b)); a=a+[Fraction(0)]*(n-len(a)); b=b+[Fraction(0)]*(n-len(b))
            return [a[i]+b[i] for i in range(n)]
        def scalar_mul(c,p): return [c*x for x in p]
        def T(p): # t p'
            return mul_t(deriv(p))
        def Q(p): # p + (t-1)p'
            # p + t p' - p'
            return add(p, add(mul_t(deriv(p)), [ -x for x in (deriv(p)+[Fraction(0)]) ])) if False else None
        # Q properly:
        def Qop(p):
            dp=deriv(p)
            tdp=mul_t(dp)
            # subtract dp (align)
            n=max(len(p),len(tdp),len(dp))
            pp=p+[Fraction(0)]*(n-len(p))
            ttdp=tdp+[Fraction(0)]*(n-len(tdp))
            dpp=dp+[Fraction(0)]*(n-len(dp))
            return [pp[i]+ttdp[i]-dpp[i] for i in range(n)]
        # recurrence
        # H_{n+1} = sum_{r=1}^{j-2} C(n,r) H_{n-r+1} + C(n,j-1)[H_{n-j+2}+ t Q(H_{n-j+1})] + 2 T H_n
        lhs=get(n+1)
        rhs=[Fraction(0)]*1
        for r in range(1, j-1):  # 1..j-2
            if r<=n:
                term=scalar_mul(math.comb(n,r), get(n-r+1))
                rhs=add(rhs,term)
        # lagged
        if n-j+2 >=0 and n-j+1>=0:
            Hlag=get(n-j+2)
            Hq=get(n-j+1)
            tQ=mul_t(Qop(Hq))
            bracket=add(Hlag, tQ)
            term=scalar_mul(math.comb(n,j-1), bracket)
            rhs=add(rhs,term)
        # 2 T H_n
        term=scalar_mul(2, T(get(n)))
        rhs=add(rhs,term)
        # pad
        N=max(len(lhs),len(rhs))
        lhs=lhs+[Fraction(0)]*(N-len(lhs))
        rhs=rhs+[Fraction(0)]*(N-len(rhs))
        if lhs!=rhs:
            errs.append((n,lhs,rhs))
    return errs

import math
for j in [2,3,4,5]:
    e=check(j,16)
    print("j",j,"errs",len(e))
    if e: print(e[0])
