import sys, math, itertools
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from fractions import Fraction

def series_inv_one_minus_A(j,N):
    A=[Fraction(0)]*(N+1)
    for r in range(1,j):
        if r<=N: A[r]=Fraction(1, math.factorial(r))
    D=[Fraction(0)]*(N+1); D[0]=Fraction(1)
    for n in range(1,N+1):
        s=Fraction(0)
        for r in range(1,min(j,n+1)):
            if r<=j-1: s+=A[r]*D[n-r]
        D[n]=s
    return D
def series_C(j,N):
    exp=[Fraction(1,math.factorial(r)) for r in range(N+1)]
    C=[Fraction(0)]*(N+1)
    for r in range(N+1):
        if r>=j: C[r]=exp[r]
    return C
def mul_series(a,b,N):
    c=[Fraction(0)]*(N+1)
    for i in range(N+1):
        for k in range(N+1-i): c[i+k]+=a[i]*b[k]
    return c

def det_frac(M):
    n=len(M); M=[row[:] for row in M]; det=Fraction(1)
    for col in range(n):
        piv=None
        for r in range(col,n):
            if M[r][col]!=0: piv=r;break
        if piv is None: return Fraction(0)
        if piv!=col: M[col],M[piv]=M[piv],M[col]; det=-det
        inv=M[col][col]; det*=inv
        for r in range(col+1,n):
            factor=M[r][col]/inv
            for c in range(col,n): M[r][c]-=factor*M[col][c]
    return det

def toeplitz(seq):
    def T(i,j):
        if i>=j and (i-j)<len(seq): return seq[i-j]
        return Fraction(0)
    return T

def minor_neg(seq, maxn, size):
    T=toeplitz(seq)
    for rows in itertools.combinations(range(maxn), size):
        for cols in itertools.combinations(range(maxn), size):
            M=[[T(r,c) for c in cols] for r in rows]
            det=det_frac(M)
            if det < 0:
                return rows, cols, det
    return None

for j in [2,3,4]:
    N=16
    C=series_C(j,N); D=series_inv_one_minus_A(j,N)
    g=mul_series(C,D,N)
    # shift: remove leading zeros -> divide by z^j
    gg=[g[r] for r in range(j, N+1)]
    print('j',j,'shifted len',len(gg))
    for size in [2,3,4]:
        r=minor_neg(gg, maxn=10, size=size)
        print('   shifted g size',size,'neg:', r)
        if r: break
