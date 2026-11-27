#!/usr/bin/env python3
import random, sympy as sp
x=sp.symbols('x')
rng=random.Random(1)
for d in range(1,7):
  for trial in range(2000):
    # q has d-1 roots, p has d roots, strictly q-interlacing p.
    qr=sorted(rng.sample(range(-30,-1),d-1)) if d>1 else []
    pr=[]
    if d==1:
      pr=[rng.randint(-30,-1)]
    else:
      for i in range(d):
        lo=(-31 if i==0 else qr[i-1]+1)
        hi=(-1 if i==d-1 else qr[i]-1)
        if lo>hi: break
        pr.append(rng.randint(lo,hi))
    if len(pr)!=d or len(set(pr))!=d: continue
    q=sp.prod(x-r for r in qr) if qr else sp.Integer(1)
    p=sp.prod(x-r for r in pr)
    n=rng.randint(1,80)
    F=sp.expand((n+1)*p+2*x*sp.diff(p,x)+n*x*(q+(x-1)*sp.diff(q,x)))
    roots=sp.nroots(F,maxsteps=200)
    if any(abs(sp.im(r))>1e-8 for r in roots):
      print('counterexample',d,n,pr,qr,sp.factor(F),roots); raise SystemExit
print('no counterexample')
