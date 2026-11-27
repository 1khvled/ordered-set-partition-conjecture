import sympy as sp
t,z,w=sp.symbols('t z w')
j=2
A=z
C=sp.exp(z)-1-z
B=1-A
D=B-t*C
F=1/D
# check H_n
def On(j,n):
    # use series
    s=sp.series(F.subs([(sp.Symbol('z'),z)]),z,0,n+1).removeO()
    return sp.expand(sp.factorial(n)*s.coeff(z,n))
# easier: use recurrence
O=[sp.Integer(1)]
for n in range(1,8):
    s=sp.Integer(0)
    for r in range(1,min(j-1,n)+1): s+=sp.binomial(n,r)*O[n-r]
    for r in range(j,n+1): s+=t*sp.binomial(n,r)*O[n-r]
    O.append(sp.expand(s))
def W(p,q): return sp.expand(sp.diff(p,t)*q-p*sp.diff(q,t))
H6 = sp.expand(t*sum(sp.binomial(6,k)*O[k] for k in range(5)))
w_direct = W(H6,O[4])
print(w_direct)
