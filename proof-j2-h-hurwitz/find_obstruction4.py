"""Obstruction with deg S=deg K=2, coupled interlacing K leftmost/S rightmost."""
import sympy as sp, itertools
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
def P(roots):
    r=1
    for x in roots: r*= (1+x*t)
    return sp.expand(r)
found=[]
for s1 in range(3,20):
 for s2 in range(1,s1):
  for b in range(s2+1,s1):
   S=P([s1,s2]); T=1+b*t
   if sp.Poly(W(S,T),t).count_roots(-sp.oo,sp.oo): continue
   for k1 in range(s1+1,25):
    for k2 in range(s2+1,s1):
     K=P([k1,k2])
     if sp.Poly(W(K,S),t).count_roots(-sp.oo,sp.oo): continue
     A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
     h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
     if sp.Poly(h0,t).count_roots(-sp.oo,sp.oo) or sp.Poly(h2,t).count_roots(-sp.oo,sp.oo): continue
     for cv in (1,2,3,5,8,13):
      H=sp.expand(h0+cv*h1+cv**2*h2)
      if sp.Poly(H,t).count_roots(-sp.oo,sp.oo):
       found.append((s1,s2,b,k1,k2,cv,sp.factor(H)))
       break
     if found: break
    if found: break
   if found: break
  if found: break
 if found: break
print('found',len(found))
for f in found[:3]: print(f[:6],f[6])
