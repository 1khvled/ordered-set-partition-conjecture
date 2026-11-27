import sys, math, itertools
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from fractions import Fraction

# g_2 coeffs
exec(open(r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-structural\pf_check.py').read().split('for j in [2,3,4,5]:')[0])
j=2
N=20
C=[Fraction(0)]*(N+1)
exp=[Fraction(1,math.factorial(r)) for r in range(N+1)]
for r in range(N+1):
    if r>=j: C[r]=exp[r]
D=series_inv_one_minus_A(j,N)
g=mul_series(C,D,N)

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

for size in [2,3,4]:
    for maxn in [8,10,12]:
        r=minor_neg(g, maxn, size)
        print('g2 size',size,'maxn',maxn,'neg:', r)
        if r: break
