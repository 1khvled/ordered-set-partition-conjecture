"""Find c where the formal model H_c fails despite coupled-type interlacings."""
import sympy as sp
t,c=sp.symbols('t c')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
S=(1+20*t)*(1+12*t); T=1+16*t; K=1+21*t
A=E(S); B=t*K; alpha=U(T); beta=E(E(K))
H=sp.expand(W(A+c*B,alpha+c*beta))
for cv in range(0,13):
    h=sp.Poly(H.subs(c,cv),t)
    print('c',cv,'roots',h.count_roots(-sp.oo,sp.oo),'H(-1/10)=',h.eval(sp.Rational(-1,10)))
