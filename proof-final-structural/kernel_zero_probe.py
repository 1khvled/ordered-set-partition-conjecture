#!/usr/bin/env python3
import sympy as sp
z=sp.symbols('z')
for j in range(2,8):
    A=sum(z**r/sp.factorial(r) for r in range(1,j))
    C=sp.series(sp.exp(z)-sum(z**r/sp.factorial(r) for r in range(j)),z,0,45).removeO()
    g=sp.series(C/(1-A),z,0,45).removeO()
    print('j',j)
    for N in range(j+1,20):
        p=sp.Poly(g,z).truncate(N) if hasattr(sp.Poly(g,z),'truncate') else sp.Poly(sum(sp.expand(g).coeff(z,r)*z**r for r in range(N)),z)
        roots=sp.nroots(p,maxsteps=200)
        nonreal=[r for r in roots if abs(sp.im(r))>1e-9]
        positive=[r for r in roots if sp.re(r)>1e-9 and abs(sp.im(r))<1e-9]
        if nonreal or positive:
            print(' first failure N',N,'nonreal',len(nonreal),'positive',len(positive),'roots',roots)
            break
    else:
        print(' all tested truncations negative-real')
