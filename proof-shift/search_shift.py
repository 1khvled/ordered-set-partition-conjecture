#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from flint import ctx
from scan import ordered_rows
from probe import add,scale,mul,deriv,roots,relation

def shift(c,s):
    out=[0]*len(c)
    for m,a in enumerate(c):
        for k in range(m+1): out[k]+=a*math.comb(m,k)*s**(m-k)
    return out

def A(p,n,s): return add(scale(p,n),scale(mul([s,1],deriv(p)),2))
def B(p,s): return mul([s,1],add(p,mul([s-1,1],deriv(p))))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--nmax',type=int,default=35); ap.add_argument('--jmax',type=int,default=8); args=ap.parse_args()
    ctx.prec=192
    for s in [0,1,2,3,4,5,8,12,20]:
        counts={}; bad=[]
        for j in range(2,args.jmax+1):
            P=[shift(c,s) for _,c in ordered_rows(j,args.nmax)]
            for n in range(1,len(P)):
                a=A(P[n-1],n,s)
                b=scale(B(P[n-j],s),math.comb(n-1,j-1)) if n>=j else [0]
                rel='none' if not any(b) else relation(roots(a),roots(b))
                counts[rel]=counts.get(rel,0)+1
                if rel=='?' and len(bad)<3: bad.append((j,n,len(a)-1,len(b)-1))
        print('shift',s,counts,'first_bad',bad)
if __name__=='__main__': main()
