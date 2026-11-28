import sympy as sp
t=sp.symbols('t')
K=sp.Function('K')(t)
def E(f): return sp.expand(f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g - f*sp.diff(g,t))
h2=W(t*K, E(E(K)))
print(sp.expand(h2))
print(sp.factor(sp.expand(h2)))
