"""W(D_p, V(xq)) bracket in delta=a-u, u, Su."""
import sympy as sp
xs,u,d,up,dp = sp.symbols('x u d up dp')
a = u+d
ap = -up + dp   # a' = u' + d' ; up = S_R (we'll set up=S_R>0 later; here up means u')
B = (3*a+xs*a**2+xs*ap)*(1+xs*u) - (2+xs*a)*(2*u+xs*u**2+xs*up)
print('B(delta,u,up,dp) =', sp.collect(sp.expand(B),[dp,d]))
