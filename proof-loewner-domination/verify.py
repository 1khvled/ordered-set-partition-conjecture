# verify.py - exact/finite audits for proof-loewner-domination
from fractions import Fraction
import numpy as np

def assoc_stirling(max_m):
    A = [[Fraction(0)]*(max_m//2+2) for _ in range(max_m+1)]
    A[0][0] = Fraction(1)
    for m in range(1, max_m+1):
        for k in range(0, m//2+1):
            val = Fraction(0)
            if k <= (m-1)//2:
                val += k * A[m-1][k]
            if k >= 1 and m >= 2 and (k-1) <= (m-2)//2:
                val += (m-1) * A[m-2][k-1]
            A[m][k] = val
    return A

def rising(a,k):
    r=Fraction(1)
    for i in range(k): r*= (a+i)
    return r

def factorial(n):
    r=1
    for i in range(2,n+1): r*=i
    return r

def O_coeffs(n,A):
    coeffs=[]
    for k in range(0,n//2+1):
        s=Fraction(0)
        for m in range(2*k,n+1):
            if k <= m//2:
                s += rising(n-m+1,k)/Fraction(factorial(m))*A[m][k]
        coeffs.append(Fraction(factorial(n))*s)
    return coeffs

def falling_factorial_coeffs(d):
    max_k=len(d)-1
    S2=[[0]*(max_k+1) for _ in range(max_k+1)]
    S2[0][0]=1
    for k in range(1,max_k+1):
        for i in range(1,k+1):
            S2[k][i]=i*S2[k-1][i]+S2[k-1][i-1]
    c=[]
    for i in range(0,max_k+1):
        s=Fraction(0)
        for k in range(i,max_k+1):
            s += Fraction(S2[k][i])*d[k]
        c.append(s)
    return c

def falling_poly(i):
    p=[Fraction(1)]
    for j in range(i):
        q=[Fraction(0)]*(len(p)+1)
        for k,pk in enumerate(p):
            q[k+1]+=pk; q[k]+=-j*pk
        p=q
    return p

def K_matrix(i,j,D):
    # returns symmetric Bezoutian kernel matrix K_{ij} in monomial basis
    fi=falling_poly(i); fj=falling_poly(j)
    def get(poly,idx): return poly[idx] if idx < len(poly) else Fraction(0)
    M=[[Fraction(0)]*D for _ in range(D)]
    maxp=max(len(fj),len(fi))
    for p in range(maxp):
        for q in range(maxp):
            Npq=get(fj,p)*get(fi,q)-get(fj,q)*get(fi,p)
            if Npq==0 or p<=q: continue
            for k in range(q,p):
                r=p-1+q-k; cc=k
                if r<D and cc<D: M[r][cc]+=Npq
    # raw is symmetric for Bezoutian; return (M+M^T)/2 for safety
    S=[[ (M[r][c]+M[c][r])/2 for c in range(D)] for r in range(D)]
    return S

def mat_add(A,B): return [[A[i][j]+B[i][j] for j in range(len(A))] for i in range(len(A))]
def mat_scale(A,s): return [[A[i][j]*s for j in range(len(A))] for i in range(len(A))]

def is_sym(A):
    D=len(A)
    for i in range(D):
        for j in range(D):
            if A[i][j]!=A[j][i]: return False
    return True

def to_float_mat(M):
    D=len(M)
    mx=max(abs(float(M[i][j])) for i in range(D) for j in range(D))
    if mx==0: mx=1.0
    A=np.array([[float(M[i][j])/mx for j in range(D)] for i in range(D)],dtype=float)
    A=(A+A.T)/2
    return A,mx

print('=== identity audits ===')
# K_{i,i+1} = f_i f_i^T
for i in range(5):
    D=i+1
    K=K_matrix(i,i+1,D)
    fi=falling_poly(i)
    vv=[[ (fi[r] if r<len(fi) else Fraction(0))*(fi[c] if c<len(fi) else Fraction(0)) for c in range(D)] for r in range(D)]
    res=max(abs(K[r][c]-vv[r][c]) for r in range(D) for c in range(D))
    print(f'K_{i},{i+1}-v_iv_i^T residual', res)

# K_{0,2}
K02=K_matrix(0,2,2)
print('K02', K02)

print('\n=== decomposition audit ===')
Aassoc=assoc_stirling(24)
for n in range(2,12):
    dP=O_coeffs(n,Aassoc); dQ=O_coeffs(n+1,Aassoc)
    cP=falling_factorial_coeffs(dP); cQ=falling_factorial_coeffs(dQ)
    D=max(len(cP)-1,len(cQ)-1)
    while len(cP)<D+1: cP.append(Fraction(0))
    while len(cQ)<D+1: cQ.append(Fraction(0))
    Delta={}
    for i in range(D+1):
        for j in range(i+1,D+1):
            Delta[(i,j)]=cP[i]*cQ[j]-cP[j]*cQ[i]
    M=[[Fraction(0)]*D for _ in range(D)]
    for i in range(D+1):
        for j in range(i+1,D+1):
            M=mat_add(M, mat_scale(K_matrix(i,j,D), Delta[(i,j)]))
    # direct Bezoutian
    # B(P,Q)(x,y)=(Q(x)P(y)-Q(y)P(x))/(x-y) in monomial basis
    # build P,Q as monomial polys from falling
    def to_mono(c):
        p=[Fraction(0)]*(len(c))
        for i,ci in enumerate(c):
            fi=falling_poly(i)
            for k,fk in enumerate(fi):
                if k < len(p): p[k]+=ci*fk
                else: p.extend([Fraction(0)]*(k-len(p)+1)); p[k]+=ci*fk
        return p
    Pm=to_mono(cP); Qm=to_mono(cQ)
    def get(p,idx): return p[idx] if idx < len(p) else Fraction(0)
    B=[[Fraction(0)]*D for _ in range(D)]
    maxp=max(len(Pm),len(Qm))
    for pidx in range(maxp):
        for qidx in range(maxp):
            Npq=get(Qm,pidx)*get(Pm,qidx)-get(Qm,qidx)*get(Pm,pidx)
            if Npq==0 or pidx<=qidx: continue
            for k in range(qidx,pidx):
                r=pidx-1+qidx-k; cc=k
                if r<D and cc<D: B[r][cc]+=Npq
    Bsym=[[ (B[r][c]+B[c][r])/2 for c in range(D)] for r in range(D)]
    res=max(abs(M[r][c]-Bsym[r][c]) for r in range(D) for c in range(D))
    print(f'n={n} D={D} decomp residual={res} symM={is_sym(M)}')

print('\n=== domination finite evidence ===')
for n in range(2,22):
    dP=O_coeffs(n,Aassoc); dQ=O_coeffs(n+1,Aassoc)
    cP=falling_factorial_coeffs(dP); cQ=falling_factorial_coeffs(dQ)
    D=max(len(cP)-1,len(cQ)-1)
    while len(cP)<D+1: cP.append(Fraction(0))
    while len(cQ)<D+1: cQ.append(Fraction(0))
    M=[[Fraction(0)]*D for _ in range(D)]
    Aadj=[[Fraction(0)]*D for _ in range(D)]
    Hlag=[[Fraction(0)]*D for _ in range(D)]
    for i in range(D+1):
        for j in range(i+1,D+1):
            K=K_matrix(i,j,D)
            w=cP[i]*cQ[j]-cP[j]*cQ[i]
            term=mat_scale(K,w)
            M=mat_add(M,term)
            if j==i+1: Aadj=mat_add(Aadj,term)
            else: Hlag=mat_add(Hlag,term)
    Mf,_=to_float_mat(M); Af,_=to_float_mat(Aadj); Hf,_=to_float_mat(Hlag)
    minM=np.linalg.eigvalsh(Mf)[0]; minA=np.linalg.eigvalsh(Af)[0]; minH=np.linalg.eigvalsh(Hf)[0]
    print(f'n={n:2d} D={D} minM={minM:.3e} minA={minA:.3e} minH={minH:.3e}')

