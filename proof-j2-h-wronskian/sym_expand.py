import sympy as sp
t,c=sp.symbols('t c')
S=sp.Function('S')(t)
K=sp.Function('K')(t)
# define operators
def E(f): return sp.expand(f + t*sp.diff(f,t))
def T_op(f): return sp.expand(sp.diff(t*sp.diff(f,t),t))
def W(f,g): return sp.expand(sp.diff(f,t)*g - f*sp.diff(g,t))

TS=T_op(S)
Pn=sp.expand(E(S)+c*t*K)
Tn=sp.expand((E(TS)+TS) + c*E(E(K)))
# verify
Pn2=sp.expand(S + t*sp.diff(S,t) + c*t*K)
Tn2=T_op(Pn2)
print(sp.expand(Tn-Tn2))

H=W(Pn2,Tn2)
# expand H in terms of derivatives of S,K
# collect by c
Hc=sp.collect(sp.expand(H),c)
print('c^2 coeff has K terms only?',)
# print coefficients
for i in range(3):
    ci=Hc.coeff(c,i)
    print('c^'+str(i)+':', sp.simplify(ci).count_ops())
