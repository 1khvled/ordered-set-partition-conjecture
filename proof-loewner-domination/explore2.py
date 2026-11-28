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
            q[k+1]+=pk
            q[k]+=-j*pk
        p=q
    return p

def K_matrix(i,j,D):
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
                if r<D and cc<D:
                    M[r][cc]+=Npq
    # symmetrize? should be symmetric
    for r in range(D):
        for c in range(r+1,D):
            if M[r][c]!=M[c][r]:
                # should not happen
                pass
    return M

def mat_scale(A,s): return [[A[i][j]*s for j in range(len(A))] for i in range(len(A))]
def mat_add(A,B): return [[A[i][j]+B[i][j] for j in range(len(A))] for i in range(len(A))]

def min_eig(M):
    A=np.array([[float(x) for x in row] for row in M],dtype=float)
    A=(A+A.T)/2
    w=np.linalg.eigvalsh(A)
    return w[0], w[-1]

max_n=16
Aassoc=assoc_stirling(max_n+2)
for n in range(2,max_n):
    dP=O_coeffs(n,Aassoc); dQ=O_coeffs(n+1,Aassoc)
    cP=falling_factorial_coeffs(dP); cQ=falling_factorial_coeffs(dQ)
    D=max(len(cP)-1,len(cQ)-1)
    while len(cP)<D+1: cP.append(Fraction(0))
    while len(cQ)<D+1: cQ.append(Fraction(0))
    Delta=[[Fraction(0)]*(D+1) for _ in range(D+1)]
    for i in range(D+1):
        for j in range(i+1,D+1):
            Delta[i][j]=cP[i]*cQ[j]-cP[j]*cQ[i]
    M=[[Fraction(0)]*D for _ in range(D)]
    Aadj=[[Fraction(0)]*D for _ in range(D)]
    Hlag=[[Fraction(0)]*D for _ in range(D)]
    for i in range(D+1):
        for j in range(i+1,D+1):
            K=K_matrix(i,j,D)
            term=mat_scale(K,Delta[i][j])
            M=mat_add(M,term)
            if j==i+1: Aadj=mat_add(Aadj,term)
            else: Hlag=mat_add(Hlag,term)
    minM,maxM=min_eig(M)
    minA,maxA=min_eig(Aadj)
    minH,maxH=min_eig(Hlag)
    # ratio / margin
    print(f"n={n} D={D} minM={minM:.3e} minA={minA:.3e} minH={minH:.3e} maxH={maxH:.3e}")
