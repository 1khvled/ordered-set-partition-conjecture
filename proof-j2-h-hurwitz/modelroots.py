import sympy as sp
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
S=(1+20*t)*(1+12*t); T=1+16*t; K=1+21*t
A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
for cv in (1,2,3):
    H=sp.Poly(sp.expand(W(A+cv*B,alpha+cv*beta)),t)
    print('c',cv,'roots',[sp.N(r) for r in sp.nroots(H)])
