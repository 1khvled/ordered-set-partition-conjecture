import sympy as sp
t=sp.symbols('t')
P=sp.Function('P')(t)
def E(f): return sp.expand(f+t*sp.diff(f,t))
def T_op(f): return sp.expand(sp.diff(t*sp.diff(f,t),t))
def W(f,g): return sp.expand(sp.diff(f,t)*g - f*sp.diff(g,t))
L=sp.expand(sp.diff(P,t)**2 - P*sp.diff(P,t,2))
H=W(P,T_op(P))
print(sp.expand(H - (E(L) - P*sp.diff(P,t,2))))
