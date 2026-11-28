"""Generic operator Wronskian identities (exact symbolic)."""
import sympy as sp
t=sp.symbols('t')
S=sp.Function('S')(t); K=sp.Function('K')(t)
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))

print('W(tK,E^2K) =', sp.factor(W(t*K,E(E(K)))))
print('W(ES,U T) =', sp.factor(W(E(S),U(sp.Function("T")(t)))))
print('W(S,T) =', sp.factor(W(S,sp.Function("T")(t))))
