import sympy as sp
z,t=sp.symbols('z t')
j=2
B=z+ t*(sp.exp(z)-1-z)
F=1/(1-B)
Fz=sp.simplify(sp.diff(F,z))
w=z**(j-1)/sp.factorial(j-1)
rhs=sp.simplify((2+(t-1)*w)*F**2 - F)
print("Fz-rhs:", sp.simplify(Fz-rhs))
