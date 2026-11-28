"""Search model with K leftmost, coupled interlacings, and H_c failing."""
import sympy as sp, itertools
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
found=[]
for a2,a1 in itertools.combinations(range(2,25),2):
  for b in range(a2+1,a1):
    S=(1+a1*t)*(1+a2*t); T=1+b*t
    if sp.Poly(W(S,T),t).count_roots(-sp.oo,sp.oo): continue
    for kv in range(1,a2):  # K root -1/kv left of -1/a2
      K=1+kv*t
      A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
      h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
      for cv in (1,2,3,5,8):
        H=sp.expand(h0+cv*h1+cv**2*h2)
        if sp.Poly(H,t).count_roots(-sp.oo,sp.oo):
          found.append((a2,a1,b,kv,cv,sp.factor(H)))
          break
      if found: break
    if found: break
  if found: break
print('found',len(found))
for f in found[:3]: print(f[:5],f[5])
