"""Search exact rational obstruction with coupled-type interlacings satisfied."""
import sympy as sp, itertools
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
def polyroots_sign(f,g):
    # check W(f,g) no real roots and sign at 0
    h=sp.Poly(W(f,g),t)
    return h.count_roots(-sp.oo,sp.oo)==0, h.eval(0)>0
found=[]
for a1,a2 in itertools.combinations(range(1,16),2):
  for b in range(a2+1,a1):
    for k in range(a1+1,20):
      S=(1+a1*t)*(1+a2*t); T=(1+b*t); K=(1+k*t)
      # S surrounds T: W(S,T)>0? determine sign via check
      okST,posST=polyroots_sign(S,T)
      if not okST: continue
      # K interlaces S appropriately: K leftmost => W(K,S) sign? just require no real roots
      okKS,_=polyroots_sign(K,S)
      if not okKS: continue
      A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
      h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
      D=sp.expand(4*h0*h2-h1**2)
      # certificate fails if h2 has real roots or D has real roots
      r2=sp.Poly(h2,t).count_roots(-sp.oo,sp.oo)
      rD=sp.Poly(D,t).count_roots(-sp.oo,sp.oo)
      if r2 or rD:
        found.append((a1,a2,b,k,r2,rD,sp.factor(h2),sp.factor(D)))
        if len(found)>=5: break
    if len(found)>=5: break
  if len(found)>=5: break
print('found',len(found))
for f in found:
  print('DATA',f[0:4],'h2 roots',f[4],'D roots',f[5])
  print(' h2=',f[6])
  print(' D=',f[7])
