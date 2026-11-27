#!/usr/bin/env python3
from __future__ import annotations
import itertools
import sympy as sp
from scan import ordered_rows

def add(dst,c,s):
    if len(dst)<len(c): dst.extend([0]*(len(c)-len(dst)))
    for i,a in enumerate(c): dst[i]+=s*a

def search(seq,label,lag=2,xpow=2,mpow=3,fit=25,validate=50):
    basis=list(itertools.product(range(lag+1),range(xpow+1),range(mpow+1)))
    eq=[]
    for m in range(lag,fit+1):
        for k in range(max(len(c) for c in seq[m-lag:m+1])+1):
            row=[]
            for l,p,a in basis:
                c=seq[m-l]
                row.append(m**a*(c[k+p] if k+p<len(c) else 0))
            eq.append(row)
    ns=sp.Matrix(eq).nullspace()
    print(label,'basis',len(basis),'nullity',len(ns))
    for v in ns[:5]:
        terms=[(basis[i],v[i]) for i in range(len(basis)) if v[i]]
        bad=[]
        for m in range(fit+1,validate+1):
            acc=[0]*(len(seq[m])+10)
            for (l,p,a),coef in terms: add(acc,seq[m-l],int(coef)*m**a)
            # shift by x^p was not applied above; redo correctly
            acc=[0]*(len(seq[m])+10)
            for (l,p,a),coef in terms:
                c=[0]*p+seq[m-l]
                add(acc,c,int(coef)*m**a)
            if any(acc): bad.append(m)
        print('terms',len(terms),'bad',bad[:3],terms)

def main():
    rows=[c for _,c in ordered_rows(2,80)]
    even=[rows[2*m] for m in range(41)]
    odd=[rows[2*m+1] for m in range(40)]
    search(even,'even')
    search(odd,'odd')
if __name__=='__main__': main()
