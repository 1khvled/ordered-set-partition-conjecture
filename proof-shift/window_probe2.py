#!/usr/bin/env python3
from __future__ import annotations
from flint import ctx
from scan import ordered_rows
from probe import add,scale,mul,deriv,roots,relation

def R(p,n): return add(scale(p,n),scale(mul([0,1],deriv(p)),2))
def main():
    ctx.prec=192
    for j in range(2,16):
        rows=[c for _,c in ordered_rows(j,100)]
        bad=[]; total=0
        for n in range(1,len(rows)):
            rpoly=R(rows[n-1],n)
            for r in range(2,j):
                if r>n: continue
                rel=relation(roots(rows[n-r]),roots(rpoly)); total+=1
                if rel=='?': bad.append((n,r,len(rows[n-r])-1,len(rpoly)-1,rel))
        print('j',j,'pairs',total,'bad',len(bad),bad[:8])
if __name__=='__main__': main()
