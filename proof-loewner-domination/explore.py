from fractions import Fraction
import math

# O_{n,2}: C(z)=e^z-1-z, F=1/(1-z-tC), O_n(t)=n![z^n]F
# d_{n,k} = [t^k]O_n(1,t) = n! * sum_{m=2k}^{n} (n-m+1 rising k / m!) * A_{m,k}
# where A_{m,k} = {m brace k}_{>=2} associated Stirling 2nd kind
# A_{m,k}=k*A_{m-1,k}+(m-1)*A_{m-2,k-1}, A_{0,0}=1

def assoc_stirling(max_m):
    # A[m][k]
    A = [[Fraction(0)]*(max_m//2+2) for _ in range(max_m+1)]
    A[0][0] = Fraction(1)
    for m in range(1, max_m+1):
        for k in range(0, m//2+1):
            val = Fraction(0)
            if k >= 0 and k <= (m-1)//2:
                val += k * A[m-1][k]
            if k >= 1 and m >= 2 and (k-1) <= (m-2)//2:
                val += (m-1) * A[m-2][k-1]
            A[m][k] = val
    return A

def rising(a, k):
    r = Fraction(1)
    for i in range(k):
        r *= (a+i)
    return r

def factorial(n):
    r = 1
    for i in range(2,n+1): r*=i
    return r

def O_coeffs(n, A):
    # d_{n,k} = n! * sum_{m=2k}^{n} rising(n-m+1,k)/m! * A_{m,k}
    coeffs = []
    for k in range(0, n//2+1):
        s = Fraction(0)
        for m in range(2*k, n+1):
            if k <= m//2:
                s += rising(n-m+1, k) / Fraction(factorial(m)) * A[m][k]
        coeffs.append(Fraction(factorial(n)) * s)
    return coeffs  # [d_0,...,d_deg]

def falling_factorial_coeffs(d):
    # O = sum_k d_k t^k = sum_i c_i (t)_i, c_i = sum_{k>=i} S2(k,i)*d_k
    # compute S2(k,i)
    max_k = len(d)-1
    S2 = [[Fraction(0)]*(max_k+1) for _ in range(max_k+1)]
    S2[0][0] = Fraction(1)
    for k in range(1, max_k+1):
        for i in range(1, k+1):
            S2[k][i] = i*S2[k-1][i] + S2[k-1][i-1]
    c = []
    for i in range(0, max_k+1):
        s = Fraction(0)
        for k in range(i, max_k+1):
            s += S2[k][i] * d[k]
        c.append(s)
    return c

def falling_poly(i):
    # (t)_i as polynomial coeffs in monomial basis [a0,...,ai]
    # (t)_0=1, (t)_{i+1}=(t-i)(t)_i
    p = [Fraction(1)]
    for j in range(i):
        # multiply by (t-j)
        q = [Fraction(0)]*(len(p)+1)
        for k,pk in enumerate(p):
            q[k+1] += pk
            q[k] += -j*pk
        p = q
    return p

def poly_add(a,b):
    n=max(len(a),len(b))
    r=[Fraction(0)]*n
    for i,x in enumerate(a): r[i]+=x
    for i,x in enumerate(b): r[i]+=x
    return r

def poly_sub(a,b):
    n=max(len(a),len(b))
    r=[Fraction(0)]*n
    for i,x in enumerate(a): r[i]+=x
    for i,x in enumerate(b): r[i]-=x
    return r

def poly_mul(a,b):
    r=[Fraction(0)]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        for j,bj in enumerate(b):
            r[i+j]+=ai*bj
    return r

def K_kernel(i,j, D):
    # K_{ij}(x,y) = (f_j(x)f_i(y)-f_j(y)f_i(x))/(x-y)
    # Return D x D coefficient matrix in monomial basis 1,x,...,x^{D-1}
    fi = falling_poly(i)
    fj = falling_poly(j)
    # numerator N(x,y)=fj(x)fi(y)-fj(y)fi(x)
    # as matrix N_{p,q} coeff of x^p y^q
    deg = max(len(fi),len(fj))-1
    size = max(len(fi)+len(fj)-2, 1)
    # build matrix of size D x D
    M = [[Fraction(0)]*D for _ in range(D)]
    # N(x,y)=sum_{p,q}(fj_p*fi_q - fj_q*fi_p) x^p y^q
    # divide by (x-y): K_{r,s} = sum_{t>=0} N_{r+1+t, s-t} - ? better use formula
    # (x^p y^q - x^q y^p)/(x-y) = sum_{k=0}^{p-q-1} x^{p-1-k} y^{q+k} for p>q
    for p in range(len(fj)):
        for q in range(len(fi)):
            coeff = fj[p]*fi[q]
            if coeff != 0:
                # add coeff * (x^p y^q - x^q y^p)/(x-y)
                add_kernel_term(M, p, q, coeff, D)
            coeff2 = -fj[q if q < len(fj) else 0]*fi[p if p < len(fi) else 0]
    # wait, easier: N_{p,q}=fj[p]*fi[q]-fj[q]*fi[p] with out of range 0
    # Let's do directly
    return M

def add_kernel_term(M, p, q, coeff, D):
    if p == q:
        return
    if p < q:
        add_kernel_term(M, q, p, -coeff, D)
        return
    # p>q: (x^p y^q - x^q y^p)/(x-y) = sum_{k=q}^{p-1} x^{p-1+q-k} y^k
    for k in range(q, p):
        r = p-1+q-k
        c = k
        if r < D and c < D:
            M[r][c] += coeff

def K_matrix(i,j,D):
    fi = falling_poly(i)
    fj = falling_poly(j)
    def get(poly, idx):
        return poly[idx] if idx < len(poly) else Fraction(0)
    M = [[Fraction(0)]*D for _ in range(D)]
    maxp = max(len(fj),len(fi))
    maxq = maxp
    for p in range(maxp):
        for q in range(maxq):
            Npq = get(fj,p)*get(fi,q) - get(fj,q)*get(fi,p)
            if Npq == 0: continue
            if p == q: continue
            if p > q:
                # add Npq * sum_{k=q}^{p-1} x^{p-1+q-k} y^k
                for k in range(q,p):
                    r = p-1+q-k
                    cc = k
                    if r < D and cc < D:
                        M[r][cc] += Npq
            # if p<q, the symmetric term (q,p) will handle it
    return M

def mat_add(A,B):
    D=len(A)
    return [[A[i][j]+B[i][j] for j in range(D)] for i in range(D)]

def mat_scale(A,s):
    D=len(A)
    return [[A[i][j]*s for j in range(D)] for i in range(D)]

def mat_sub(A,B):
    D=len(A)
    return [[A[i][j]-B[i][j] for j in range(D)] for i in range(D)]

def is_symmetric(M):
    D=len(M)
    for i in range(D):
        for j in range(D):
            if M[i][j]!=M[j][i]: return False
    return True

def leading_principal_minors(M):
    # return list of det of top-left kxk
    D=len(M)
    res=[]
    for k in range(1,D+1):
        # compute det by bareiss or fraction-free? use Fraction gaussian
        A=[[M[i][j] for j in range(k)] for i in range(k)]
        res.append(det(A))
    return res

def det(A):
    n=len(A)
    M=[row[:] for row in A]
    sign=Fraction(1)
    for col in range(n):
        pivot=None
        for row in range(col,n):
            if M[row][col]!=0:
                pivot=row; break
        if pivot is None:
            return Fraction(0)
        if pivot!=col:
            M[col],M[pivot]=M[pivot],M[col]
            sign=-sign
        pv=M[col][col]
        for row in range(col+1,n):
            if M[row][col]==0: continue
            f=M[row][col]/pv
            for j in range(col,n):
                M[row][j]-=f*M[col][j]
    d=sign
    for i in range(n): d*=M[i][i]
    return d

def eigenvalues_sym_exact_small(M):
    # only for 2x2
    if len(M)==2:
        a,b,c,d=M[0][0],M[0][1],M[1][0],M[1][1]
        # a,d + (b+c)/2? symmetric
        tr=a+d
        disc=(a-d)*(a-d)+4*b*c
        return (tr, disc)
    return None

def mat_is_psd(M):
    # check all principal minors? or leading principal minors for PD; for PSD need all principal minors
    # use Sylvester's criterion for PD; for PSD check eigenvalues via numpy? No, use exact: check all principal minors >=0
    D=len(M)
    # check all principal minors
    for mask in range(1, 1<<D):
        idx=[i for i in range(D) if mask&(1<<i)]
        sub=[[M[i][j] for j in idx] for i in idx]
        if det(sub) < 0:
            return False
    return True

def mat_is_pd(M):
    minors=leading_principal_minors(M)
    return all(m>0 for m in minors)

def min_eigenvalue_numeric(M):
    import numpy as np
    A=np.array([[float(x) for x in row] for row in M], dtype=float)
    w=np.linalg.eigvalsh(A)
    return float(w[0]), float(w[-1])

max_n=10
Aassoc = assoc_stirling(max_n+2)
for n in range(2, max_n):
    dP = O_coeffs(n, Aassoc)
    dQ = O_coeffs(n+1, Aassoc)
    cP = falling_factorial_coeffs(dP)
    cQ = falling_factorial_coeffs(dQ)
    D = max(len(cP)-1, len(cQ)-1)
    # pad
    while len(cP)<D+1: cP.append(Fraction(0))
    while len(cQ)<D+1: cQ.append(Fraction(0))
    # Delta
    Delta = [[Fraction(0)]*(D+1) for _ in range(D+1)]
    for i in range(D+1):
        for j in range(i+1,D+1):
            Delta[i][j]=cP[i]*cQ[j]-cP[j]*cQ[i]
    # M
    M = [[Fraction(0)]*D for _ in range(D)]
    Aadj = [[Fraction(0)]*D for _ in range(D)]
    Hlag = [[Fraction(0)]*D for _ in range(D)]
    for i in range(D+1):
        for j in range(i+1,D+1):
            K = K_matrix(i,j,D)
            term = mat_scale(K, Delta[i][j])
            M = mat_add(M, term)
            if j==i+1:
                Aadj = mat_add(Aadj, term)
            else:
                Hlag = mat_add(Hlag, term)
    sym = is_symmetric(M)
    pd = mat_is_pd(M)
    # check domination Aadj + Hlag = M; check if Aadj dominates -Hlag in some sense?
    # Loewner domination: M = Aadj + Hlag >= ? . The question says "domination of indefinite kernels by adjacent positive kernels"
    # so Aadj + Hlag should be PSD, i.e. Aadj >= -Hlag in Loewner? -Hlag may be indefinite too.
    # Direct check: M PD?
    # Also check Aadj - |Hlag| not right. Check eigenvalues of M, Aadj, Hlag
    try:
        minM, maxM = min_eigenvalue_numeric(M)
        minA, maxA = min_eigenvalue_numeric(Aadj)
        minH, maxH = min_eigenvalue_numeric(Hlag)
    except Exception as e:
        minM=maxM=minA=maxA=minH=maxH=float('nan')
    # check Aadj + t*Hlag? no. Check if M >= 0.
    print(f"n={n} D={D} sym={sym} M_PD={pd} minEigM={minM:.6g} minEigA={minA:.6g} minEigH={minH:.6g}")
    # print Deltas
    dels=[]
    for i in range(D+1):
        for j in range(i+1,D+1):
            if Delta[i][j]!=0:
                dels.append((i,j,Delta[i][j]))
    print("  Deltas:", dels[:6], "..." if len(dels)>6 else "")
    # check if any Delta negative
    neg=[x for x in dels if x[2]<0]
    if neg: print("  NEG DELTA", neg)
