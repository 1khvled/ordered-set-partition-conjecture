"""Broad search for a counterexample to H>0 (fast integer coefficients + scaled roots)."""
import numpy as np, math

def triangle(nmax):
    w=nmax//2+3
    A=[[0]*w for _ in range(nmax+1)]
    A[0][0]=1
    for m in range(nmax):
        for r in range(w-1):
            A[m+1][r]+=r*A[m][r]
            if r: A[m+1][r]+=m*A[m-1][r-1]
    return A

def rising(a,r):
    v=1
    for i in range(r): v*=a+i
    return v

def coeffs_P(A,n,a):
    # p[k] = (a)_{k+1} A[n][k+1]
    return [rising(a,k+1)*A[n][k+1] for k in range(len(A[n])-1)]

def coeffs_T(p):
    # T = sum k^2 p[k] t^{k-1}
    return [(k+1)**2*p[k+1] for k in range(len(p)-1)]

def W_coeffs(p,q):
    # W = p' q - p q'; return ascending coeffs
    dp=[k*p[k] for k in range(1,len(p))]
    dq=[k*q[k] for k in range(1,len(q))]
    res=[0]*(len(dp)+len(q)-1+1)
    for i,ci in enumerate(dp):
        for j,cj in enumerate(q): res[i+j]+=ci*cj
    for i,ci in enumerate(p):
        for j,cj in enumerate(dq): res[i+j]-=ci*cj
    while res and res[-1]==0: res.pop()
    return res

N=130
A=triangle(N)
worst=[]
for n in range(4,N+1):
    for a in range(1,31):
        p=coeffs_P(A,n,a)
        q=coeffs_T(p)
        h=W_coeffs(p,q)
        if len(h)<2: continue
        scale=max(abs(x) for x in h)
        sh=max(0,scale.bit_length()-900)
        hs=[x>>sh for x in h]
        sc=max(abs(x) for x in hs)
        cs=np.array([float(x)/sc for x in hs])
        roots=np.roots(cs[::-1])
        mx=max(z.real for z in roots)
        if mx>1e-8 or min(abs(z.imag) for z in roots)<1e-12:
            worst.append((n,a,mx,min(abs(z.imag) for z in roots)))
    if n%20==0: print('done n',n,flush=True)
print('suspects',len(worst))
for x in worst[:20]: print(x)
