"""Bracket for W(P,Q) with R=xq: P=p+V p, Q=V R."""
import sympy as sp
xs,a,u,ap,up=sp.symbols('x a u ap up')
B = (3*a+xs*a**2+xs*ap)*(1+xs*u) - (2+xs*a)*(2*u+xs*u**2+xs*up)
print('B =', sp.factor(sp.expand(B)))
Sp,Sq=sp.symbols('Sp Sq', nonnegative=True)
print('B with ap=-Sp,up=-Sq:', sp.collect(sp.expand(B.subs({ap:-Sp,up:-Sq})),[Sp,Sq]))
