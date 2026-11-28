"""Correct exact finite search for a positive zero of W(R,(xR')')."""
import sympy as s
x=s.symbols('x'); N=60
A=[[0]*(N//2+3) for _ in range(N+1)]; A[0][0]=1
for m in range(N):
 for r in range(len(A[0])-1):
  A[m+1][r]+=r*A[m][r]
  if r: A[m+1][r]+=m*A[m-1][r-1]
def rising(a,r):
 v=1
 for i in range(r): v*=a+i
 return v
def H(n,a):
 R=s.expand(sum(rising(a,r)*A[n][r]*(-x)**(r-1) for r in range(1,len(A[0])) if A[n][r]))
 C=s.expand(s.diff(x*s.diff(R,x),x))
 return s.Poly(s.expand(s.diff(R,x)*C-R*s.diff(C,x)),x)
for a in [1,2,3,5,10]:
 for n in range(6,N+1):
  h=H(n,a); nr=h.count_roots(0,s.oo)
  if nr:
   print('CANDIDATE',n,a,'H positive roots',nr,flush=True); raise SystemExit
 print('a',a,'clear through',N,flush=True)
print('no exact candidate')
