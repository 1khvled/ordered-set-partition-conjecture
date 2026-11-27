import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')

# Compute power series coefficients (ordinary) of g_j = C_j/(1-A_j) and D_j = 1/(1-A_j)
# A_j = sum_{r=1}^{j-1} z^r/r!, C_j = e^z - sum_{r=0}^{j-1} z^r/r!, D=1/(1-A)
# Work with rational coefficients via fractions.

from fractions import Fraction

def series_exp(N):
    # coefficients of e^z as ordinary power series: 1/r!
    c=[Fraction(1, math.factorial(r)) for r in range(N+1)]
    return c

def series_inv_one_minus_A(j, N):
    # A = sum_{r=1}^{j-1} z^r/r!
    # D = 1/(1-A): solve D*(1-A)=1
    A=[Fraction(0)]*(N+1)
    for r in range(1, j):
        if r<=N: A[r]=Fraction(1, math.factorial(r))
    D=[Fraction(0)]*(N+1)
    D[0]=Fraction(1)
    for n in range(1, N+1):
        # D_n = sum_{r=1}^{min(j-1,n)} A_r * D_{n-r}  (since (1-A)D=1 => D = 1 + A D)
        s=Fraction(0)
        for r in range(1, min(j, n+1)):
            if r<=j-1:
                s += A[r]*D[n-r]
        D[n]=s
    return D

def series_C(j, N):
    exp=series_exp(N)
    C=[Fraction(0)]*(N+1)
    for r in range(N+1):
        if r < j:
            C[r]=Fraction(0)  # subtract terms r=0..j-1
        else:
            C[r]=exp[r]
    # Actually C = e^z - sum_{r=0}^{j-1} z^r/r!, so coeff at r<j is 0 for r>=1? term r=0: 1-1=0. yes.
    return C

def mul_series(a,b,N):
    c=[Fraction(0)]*(N+1)
    for i in range(N+1):
        for k in range(N+1-i):
            c[i+k]+=a[i]*b[k]
    return c

def div_series(num, den, N):
    # den[0]!=0
    res=[Fraction(0)]*(N+1)
    for n in range(N+1):
        s=num[n]
        for k in range(1, n+1):
            s -= den[k]*res[n-k]
        res[n]=s/den[0]
    return res

def coeffs_g(j,N):
    C=series_C(j,N)
    D=series_inv_one_minus_A(j,N)
    # g = C * D = C/(1-A)
    return mul_series(C,D,N)

def coeffs_D(j,N):
    return series_inv_one_minus_A(j,N)

def toeplitz_minor_nonnegative(seq, size):
    # check all minors of given size? just first few leading
    import itertools
    N=len(seq)
    # Build Toeplitz matrix T_{i,j}=seq[i-j] (i>=j) maybe; for PF sequence need all minors of T (infinite) nonneg.
    # Test leading principal minors up to size
    for sz in range(2, size+1):
        # leading sz x sz
        M=[[ (seq[i-j] if i>=j else 0) for j in range(sz)] for i in range(sz)]
        # determinant
        from fractions import Fraction
        # compute det
        import copy
        M=[row[:] for row in M]
        det=Fraction(1)
        n=sz
        for col in range(n):
            # find pivot
            piv=None
            for r in range(col,n):
                if M[r][col]!=0:
                    piv=r; break
            if piv is None:
                det=Fraction(0); break
            if piv!=col:
                M[col],M[piv]=M[piv],M[col]; det=-det
            inv=M[col][col]
            det*=inv
            for r in range(col+1,n):
                factor=M[r][col]/inv
                for c in range(col,n):
                    M[r][c]-=factor*M[col][c]
        # det already product of pivots approx; but elimination changes? We did LU without dividing rows, det = product pivots * sign. Actually after elimination, determinant = product of original pivots (with sign). We tracked det*=inv before elimination, so det = product pivots. Good.
        if det < 0:
            return sz, det
    return None

for j in [2,3,4,5]:
    N=12
    g=coeffs_g(j,N)
    D=coeffs_D(j,N)
    print('j',j)
    print('  g coeffs:', [str(x) for x in g[:8]])
    print('  D coeffs:', [str(x) for x in D[:8]])
    r=toeplitz_minor_nonnegative(g, 6)
    print('  g PF leading minor failure:', r)
    r2=toeplitz_minor_nonnegative(D, 6)
    print('  D PF leading minor failure:', r2)
