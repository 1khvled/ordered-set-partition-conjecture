#!/usr/bin/env python3
from __future__ import annotations
import math
from flint import ctx
from scan import ordered_rows
from probe import add,scale,mul,deriv,roots,relation

def A(p,n): return add(scale(p,n),scale(mul([0,1],deriv(p)),2))
def small_sum(rows,n,j):
    out=[0]
    for r in range(2,j):
        if r<=n: out=add(out,scale(rows[n-r],math.comb(n,r)))
    return out

def main():
    ctx.prec=192
    for j in range(2,16):
        rows=[c for _,c in ordered_rows(j,100)]
        counts={}; bad=[]
        for n in range(1,len(rows)):
            q=rows[n-j] if n>=j else [1]
            a=A(rows[n-1],n)
            rel=relation(roots(q),roots(a))
            counts[rel]=counts.get(rel,0)+1
            if rel=='?' and len(bad)<5: bad.append((n,len(q)-1,len(a)-1))
            if j>2:
                aa=add(a,small_sum(rows,n,j))
                rel2=relation(roots(q),roots(aa))
                counts['small_'+rel2]=counts.get('small_'+rel2,0)+1
                if rel2=='?' and len(bad)<5: bad.append(('small',n,len(q)-1,len(aa)-1))
        print('j',j,counts,'bad',bad)
if __name__=='__main__': main()
