import sympy as sp
z,t=sp.symbols('z t')
def E(j): return sum(z**r/sp.factorial(r) for r in range(j))
for j in [2,3,4]:
    w=z**(j-1)/sp.factorial(j-1)
    q2=2+(t-1)*w
    integral=sp.integrate(q2*sp.exp(-z),(z,0,z))
    invF=sp.simplify(sp.exp(z)*(1-integral))
    B=sum(c for c in [])  # placeholder
    # original B_j
    Aj=sum(z**r/sp.factorial(r) for r in range(1,j))
    Cj=sp.exp(z)-E(j)
    Bj=Aj + t*Cj
    orig=sp.simplify(1/(1/(1-Bj)))  # 1-B? just compute 1-Bj
    orig_den=sp.simplify(1 - Bj)
    print(f"j={j}")
    print("  new 1/F:", sp.simplify(invF))
    print("  orig 1-B:", sp.simplify(orig_den))
    print("  diff:", sp.simplify(invF-orig_den))
PY