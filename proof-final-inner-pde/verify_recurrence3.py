import sys, math
sys.path.insert(0, r"C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-inner")
from explore_H import inner_array
from fractions import Fraction

def H_coeffs(j,n,nmax=20):
    h = inner_array(j, nmax)
    row = h[n]
    deg = n//j
    fact = math.factorial(n)
    return [fact*row[m] for m in range(deg+1)]

def check(j,nmax):
    def get(N): return H_coeffs(j,N,nmax)
    def deriv(p):
        return [ (k+1)*p[k+1] for k in range(len(p)-1) ] if len(p)>1 else [Fraction(0)]
    def mul_t(p): return [Fraction(0)]+p
    def add(a,b):
        n=max(len(a),len(b)); a=a+[Fraction(0)]*(n-len(a)); b=b+[Fraction(0)]*(n-len(b))
        return [a[i]+b[i] for i in range(n)]
    def scalar_mul(c,p): return [c*x for x in p]
    def T(p): return mul_t(deriv(p))
    def Qop(p):
        dp=deriv(p); tdp=mul_t(dp)
        n=max(len(p),len(tdp),len(dp))
        pp=p+[Fraction(0)]*(n-len(p)); ttdp=tdp+[Fraction(0)]*(n-len(tdp)); dpp=dp+[Fraction(0)]*(n-len(dp))
        return [pp[i]+ttdp[i]-dpp[i] for i in range(n)]
    errs=[]
    for n in range(1, nmax-1):
        lhs=get(n+1)
        rhs=[Fraction(0)]
        for r in range(1, j-1):  # 1..j-2
            if r<=n:
                term=scalar_mul(math.comb(n,r), get(n-r+1))
                rhs=add(rhs,term)
        if n-j+2 >=0 and n-j+1>=0:
            Hlag=get(n-j+2); Hq=get(n-j+1)
            tQ=mul_t(Qop(Hq))
            bracket=add(Hlag, tQ)
            term=scalar_mul(math.comb(n,j-1), bracket)
            rhs=add(rhs,term)
        term=scalar_mul(2, T(get(n)))
        rhs=add(rhs,term)
        N=max(len(lhs),len(rhs))
        lhs=lhs+[Fraction(0)]*(N-len(lhs)); rhs=rhs+[Fraction(0)]*(N-len(rhs))
        if lhs!=rhs:
            errs.append((n))
    return errs

for j in [2,3,4,5]:
    e=check(j,16)
    print("j",j,"errs",len(e))
    if e: print(e[:5])
