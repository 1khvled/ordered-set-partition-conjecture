from fractions import Fraction
import numpy as np

# compute M, Aadj, Hlag explicitly for n=5 (D=3)
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
    # make symmetric
    Ms=[[M[i][j] for j in range(D)] for i in range(D)]
    for r in range(D):
        for c in range(r+1,D):
            v=Ms[r][c]+Ms[c][r]
            Ms[r][c]=v; Ms[c][r]=v
    return Ms

Aassoc=assoc_stirling(15)
for n in [5]:
    dP=O_coeffs(n,Aassoc); dQ=O_coeffs(n+1,Aassoc)
    cP=falling_factorial_coeffs(dP); cQ=falling_factorial_coeffs(dQ)
    D=max(len(cP)-1,len(cQ)-1)
    while len(cP)<D+1: cP.append(Fraction(0))
    while len(cQ)<D+1: cQ.append(Fraction(0))
    print('cP',cP)
    print('cQ',cQ)
    Delta={}
    for i in range(D+1):
        for j in range(i+1,D+1):
            Delta[(i,j)]=cP[i]*cQ[j]-cP[j]*cQ[i]
            print(f'Delta_{i}{j}={Delta[(i,j)]}')
    for i in range(D+1):
        for j in range(i+1,D+1):
            K=K_matrix(i,j,D)
            print(f'K_{i}{j}=')
            for row in K: print(row)
    M=[[Fraction(0)]*D for _ in range(D)]
    for i in range(D+1):
        for j in range(i+1,D+1):
            K=K_matrix(i,j,D)
            for r in range(D):
                for c in range(D):
                    M[r][c]+=Delta[(i,j)]*K[r][c]
    print('M=')
    for row in M: print(row)
