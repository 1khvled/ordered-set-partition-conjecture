import sys, math, itertools
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from fractions import Fraction

def coeff_C(j,N):
    exp=[Fraction(1,math.factorial(r)) for r in range(N+1)]
    C=[Fraction(0)]*(N+1)
    for r in range(N+1):
        if r>=j:
            C[r]=exp[r]
    return C

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

def minor_neg(seq, maxn=9, size=2):
    T=toeplitz(seq)
    for rows in itertools.combinations(range(maxn), size):
        for cols in itertools.combinations(range(maxn), size):
            M=[[T(r,c) for c in cols] for r in rows]
            det=det_frac(M)
            if det < 0:
                return rows, cols, det
    return None

for j in [2,3,4,5]:
    C=coeff_C(j,14)
    print('j',j,'C minor neg size2:', minor_neg(C,9,2), 'size3:', minor_neg(C,9,3))
