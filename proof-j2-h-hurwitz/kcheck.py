import sympy as sp
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
S=(1+20*t)*(1+12*t); T=1+16*t
for kv in (4,6,8,10,11,13,21):
    K=1+kv*t
    A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
    h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
    D=sp.expand(4*h0*h2-h1**2)
    print('k',kv,'h2roots',sp.Poly(h2,t).count_roots(-sp.oo,sp.oo),'Droots',sp.Poly(D,t).count_roots(-sp.oo,sp.oo),'D(0)',sp.Poly(D,t).eval(0))
