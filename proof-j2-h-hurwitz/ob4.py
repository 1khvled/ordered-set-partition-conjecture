import sympy as sp
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
S=(1+20*t)*(1+12*t); T=1+16*t; K=1+4*t
A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
D=sp.expand(4*h0*h2-h1**2)
print('h0',sp.factor(h0)); print('h2',sp.factor(h2)); print('D',sp.factor(D))
for cv in (1,2,3,5):
    H=sp.expand(h0+cv*h1+cv**2*h2)
    print('c',cv,'realroots',sp.Poly(H,t).count_roots(-sp.oo,sp.oo),'roots', [sp.N(r) for r in sp.nroots(sp.Poly(H,t))])
