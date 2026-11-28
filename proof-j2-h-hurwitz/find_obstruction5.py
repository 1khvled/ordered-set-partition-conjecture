"""Correct-orientation obstruction search: roots -sA<-sB, K roots -kA<-kB with kA>sA>kB>sB."""
import sympy as sp, itertools
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
def P(rs):
    r=1
    for x in rs: r*=(1+t/x)
    return sp.cancel(r)
found=[]
for sB in range(1,10):
 for sA in range(sB+1,20):
  for b in range(sB+1,sA):
   S=P([sA,sB]); T=1+t/sp.Integer(b)
   if sp.Poly(W(S,T),t).count_roots(-sp.oo,sp.oo): continue
   for kB in range(sB+1,sA):
    for kA in range(sA+1,25):
     K=P([kA,kB])
     if sp.Poly(W(K,S),t).count_roots(-sp.oo,sp.oo): continue
     A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
     h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
     if sp.Poly(h0,t).count_roots(-sp.oo,sp.oo) or sp.Poly(h2,t).count_roots(-sp.oo,sp.oo): continue
     for cv in (1,2,3,5,8,13,21):
      H=sp.expand(h0+cv*h1+cv**2*h2)
      if sp.Poly(H,t).count_roots(-sp.oo,sp.oo):
       found.append((sA,sB,b,kA,kB,cv,sp.factor(H)))
       break
     if found: break
    if found: break
   if found: break
  if found: break
 if found: break
print('found',len(found))
for f in found[:3]: print(f[:6],f[6])
