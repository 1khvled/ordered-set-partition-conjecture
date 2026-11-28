"""Debug specific candidate obstruction."""
import sympy as sp
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
S=(1+20*t)*(1+12*t); T=1+16*t; K=1+21*t
print('W(S,T) roots',sp.Poly(W(S,T),t).count_roots(-sp.oo,sp.oo),sp.factor(W(S,T)))
A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
D=sp.expand(4*h0*h2-h1**2)
print('h0',sp.factor(h0),'roots',sp.Poly(h0,t).count_roots(-sp.oo,sp.oo))
print('h2',sp.factor(h2),'roots',sp.Poly(h2,t).count_roots(-sp.oo,sp.oo))
print('D roots',sp.Poly(D,t).count_roots(-sp.oo,sp.oo),'D',sp.factor(D))
