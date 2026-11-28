import sympy as sp
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
S=(1+18*t)*(1+4*t); T=1+5*t; K=1+t
A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
H=sp.expand(h0+5*h1+25*h2)
for name,p in [('W(S,T)',W(S,T)),('W(K,S)',W(K,S)),('h0',h0),('h2',h2),('H',H)]:
    print(name,'=',sp.factor(p),'realroots',sp.Poly(p,t).count_roots(-sp.oo,sp.oo),'at0',sp.Poly(p,t).eval(0))
print('H roots',[sp.N(r) for r in sp.nroots(sp.Poly(H,t))])
