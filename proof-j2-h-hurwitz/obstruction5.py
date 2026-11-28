"""Exact obstruction with the correct K-leftmost coupled orientation."""
import sympy as sp
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
q=sp.Rational
S=(1+t/q(8))*(1+t); T=1+t/q(7); K=(1+t/q(9))*(1+t/q(3))
A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
D=sp.expand(4*h0*h2-h1**2)
H=sp.expand(h0+h1+h2)
for name,p in [('W(S,T)',W(S,T)),('W(K,S)',W(K,S)),('h0',h0),('h2',h2),('D',D),('H(c=1)',H)]:
    print(name,'realroots',sp.Poly(p,t).count_roots(-sp.oo,sp.oo),'at0',sp.Poly(p,t).eval(0))
print('H roots',[sp.N(r) for r in sp.nroots(sp.Poly(H,t))])
print('H=',sp.factor(H))
