import sympy as sp
t=sp.symbols('t')
S6=502+830*t+90*t**2
O4=24+45*t+6*t**2
w=sp.expand(sp.diff(S6,t)*O4 - S6*sp.diff(O4,t))
print(w)
print(sp.factor(w))
