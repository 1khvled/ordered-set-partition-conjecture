"""Relaxed search: S deg2 surrounds T deg1, K deg1 arbitrary."""
import sympy as sp, itertools
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
found=[]
for a2,a1 in itertools.combinations(range(1,30),2):
  for b in range(a2+1,a1):
    S=(1+a1*t)*(1+a2*t); T=(1+b*t)
    if sp.Poly(W(S,T),t).count_roots(-sp.oo,sp.oo): continue
    for k in range(1,25):
      if k in (a1,a2): continue
      K=(1+k*t)
      A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
      h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
      D=sp.expand(4*h0*h2-h1**2)
      r2=sp.Poly(h2,t).count_roots(-sp.oo,sp.oo)
      rD=sp.Poly(D,t).count_roots(-sp.oo,sp.oo)
      if r2 or rD:
        found.append((a1,a2,b,k,r2,rD,sp.factor(h2),sp.factor(D)))
print('found',len(found))
for f in found[:6]:
  print('DATA S=(1+%dt)(1+%dt) T=(1+%dt) K=(1+%dt): h2roots %d Droots %d'%(f[0],f[1],f[2],f[3],f[4],f[5]))
  print(' h2=',f[6]); print(' D=',f[7])
