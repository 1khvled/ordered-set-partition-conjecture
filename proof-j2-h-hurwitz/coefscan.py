"""Exact coefficient-positivity scan, broad range."""
from bigsearch import triangle, coeffs_P, coeffs_T, W_coeffs
N=140
A=triangle(N)
bad=[]
for n in range(4,N+1):
    for a in range(1,41):
        p=coeffs_P(A,n,a); q=coeffs_T(p); h=W_coeffs(p,q)
        if any(c<=0 for c in h): bad.append((n,a,min(h)))
print('bad coefficient cases',len(bad),bad[:10])
