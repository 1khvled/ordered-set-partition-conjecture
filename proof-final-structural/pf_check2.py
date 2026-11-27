import sys, math, itertools
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from fractions import Fraction

exec(open(r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-structural\pf_check.py').read().split('for j in [2,3,4,5]:')[0])

def toeplitz(seq):
    # infinite lower-triangular Toeplitz T_{i,j}=seq[i-j] if i>=j else 0
    def T(i,j):
        if i>=j and (i-j)<len(seq):
            return seq[i-j]
        return Fraction(0)
    return T

def minor_neg(seq, maxn=8, size=2):
    T=toeplitz(seq)
    for rows in itertools.combinations(range(maxn), size):
        for cols in itertools.combinations(range(maxn), size):
            # build
            M=[[T(r,c) for c in cols] for r in rows]
            # det
            det = det_frac(M)
            if det < 0:
                return rows, cols, det
    return None

def det_frac(M):
    n=len(M)
    M=[row[:] for row in M]
    det=Fraction(1)
    for col in range(n):
        piv=None
        for r in range(col,n):
            if M[r][col]!=0:
                piv=r;break
        if piv is None:
            return Fraction(0)
        if piv!=col:
            M[col],M[piv]=M[piv],M[col]; det=-det
        inv=M[col][col]
        det*=inv
        for r in range(col+1,n):
            factor=M[r][col]/inv
            for c in range(col,n):
                M[r][c]-=factor*M[col][c]
    return det

for j in [2,3,4]:
    N=14
    g=coeffs_g(j,N)
    D=coeffs_D(j,N)
    print('j',j)
    print('  g minor neg:', minor_neg(g, maxn=9, size=2))
    print('  g minor neg size3:', minor_neg(g, maxn=9, size=3))
    print('  D minor neg:', minor_neg(D, maxn=9, size=2))
    print('  D minor neg size3:', minor_neg(D, maxn=9, size=3))
