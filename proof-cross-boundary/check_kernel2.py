import sympy as sp
t,z,w=sp.symbols('t z w')
j=2
A=z
Cz=sp.exp(z)-1-z
Bz=1-A
Dz=Bz-t*Cz
Cw=Cz.subs(z,w)
Bw=Bz.subs(z,w)
Dw=Dw_=Bw-t*Cw
# kernel
K = Cz*(Bz*Dw - t*Cw*Dz)/(Dz**2 * Dw**2)
# expand
n=6; d=2
serz = sp.series(K, z, 0, n+1).removeO()
coeffz = serz.coeff(z,n)
serw = sp.series(coeffz, w, 0, n-d+1).removeO()
coeffw = serw.coeff(w, n-d)
res = sp.expand(sp.factorial(n)*sp.factorial(n-d)*coeffw)
print(res)
print(res == sp.expand(540*t**4 + 8100*t**3 + 40728*t**2 + 39840*t + 12408))
