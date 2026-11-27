import sympy as sp
z,t=sp.symbols('z t')
j=2
w=z**(j-1)/sp.factorial(j-1)
q2=2+(t-1)*w
# v = C q2 e^{-z}; choose C=-1: u = 1 + integral? Let's just solve: u' = -q2*exp(-z) (C=-1), u(0)=1
integrand=sp.integrate(-q2*sp.exp(-z),(z,0,z))
u=1+integrand
F_new=sp.simplify(sp.exp(-z)/u)
B=z+ t*(sp.exp(z)-1-z)
F_orig=1/(1-B)
print("series new:", sp.series(F_new,z,0,4))
print("series orig:", sp.series(F_orig,z,0,4))
PY