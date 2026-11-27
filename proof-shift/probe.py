#!/usr/bin/env python3
"""Probe the x=t-1 transformed universal recurrence for compatibility."""
from __future__ import annotations
import json
from pathlib import Path
from flint import ctx, fmpz_poly
from scan import ordered_rows, certified_complex_roots


def add(a,b):
    c=[0]*max(len(a),len(b))
    for i,x in enumerate(a): c[i]+=x
    for i,x in enumerate(b): c[i]+=x
    return c

def scale(a,k): return [k*x for x in a]
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def deriv(a): return [i*a[i] for i in range(1,len(a))] or [0]
def shift(c): # p(x+1)
    out=[0]*len(c)
    for m,a in enumerate(c):
        for k in range(m+1): out[k]+=a*__import__('math').comb(m,k)
    return out
def A(p,n): return add(scale(p,n),scale(mul([1,1],deriv(p)),2))
def B(p): return mul([1,1],add(p,mul([0,1],deriv(p))))

def leq(a,b): return bool(a.upper()<=b.lower())
def roots(c):
    rr=certified_complex_roots(c,192)
    if any(not r.imag.is_zero() for r,_ in rr): raise ValueError('nonreal')
    v=[r.real for r,m in rr for _ in range(m)]
    v.sort(key=lambda x:float(x.mid()))
    if any(not leq(x,y) for x,y in zip(v,v[1:])): raise ArithmeticError('overlap')
    return v

def relation(a,b):
    # Return L if a alternates left of b, R if b alternates left of a,
    # I if lower degree interlaces higher, or ?.
    d,e=len(a),len(b)
    if d==e:
        left=all(leq(a[i],b[i]) for i in range(d)) and all(leq(b[i],a[i+1]) for i in range(d-1))
        right=all(leq(b[i],a[i]) for i in range(d)) and all(leq(a[i],b[i+1]) for i in range(d-1))
        if left: return 'L'
        if right: return 'R'
        return '?'
    if d+1==e:
        return 'I' if all(leq(b[i],a[i]) and leq(a[i],b[i+1]) for i in range(d)) else '?'
    if e+1==d:
        return 'I' if all(leq(a[i],b[i]) and leq(b[i],a[i+1]) for i in range(e)) else '?'
    return 'jump'


def main():
    ctx.prec=192
    out=Path('shift-results.jsonl')
    counts={}
    examples=[]
    with out.open('w',encoding='utf-8') as f:
        for j in range(2,21):
            P=[shift(c) for _,c in ordered_rows(j,120)]
            for n in range(1,len(P)):
                a=A(P[n-1],n)
                b=scale(B(P[n-j]),__import__('math').comb(n-1,j-1)) if n>=j else [0]
                if not any(b):
                    rel='none'; ok=True
                else:
                    rel=relation(roots(a),roots(b)); ok=rel!='?'
                counts[rel]=counts.get(rel,0)+1
                rec={'j':j,'n':n,'relation':rel,'ok':ok,'deg_a':len(a)-1,'deg_b':len(b)-1}
                f.write(json.dumps(rec)+'\n')
                if not ok and len(examples)<20: examples.append(rec)
    print(counts)
    print(json.dumps(examples,indent=2))

if __name__=='__main__': main()
