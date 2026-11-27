import sympy as sp
from sympy import symbols,factorial,binomial,exp,Rational,Poly,series

z,t,u=symbols('z t u')
for j in [2,3,4,5]:
    A=sum(z**r/factorial(r) for r in range(1,j))
    C=exp(z)-sum(z**r/factorial(r) for r in range(0,j))
    # prompt's claimed identity
    lhs = 1-A-t*C
    rhs = (1+t)*(1-A - (t/(1+t))*exp(z))
    diff = sp.simplify(lhs-rhs)
    print('j',j,'prompt diff', sp.simplify(diff))
    # correct identity
    lhs2 = sp.simplify((1-u)*(1-A - (u/(1-u))*C))
    rhs2 = sp.simplify(1 - A*(1-2*u) - u*exp(z))
    print('  correct diff', sp.simplify(lhs2-rhs2))
