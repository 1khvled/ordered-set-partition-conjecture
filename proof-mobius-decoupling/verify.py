#!/usr/bin/env python3
"""Exact checks for Mobius decoupling attempt."""
from __future__ import annotations
import sympy as sp
from sympy import symbols, factorial, binomial, exp, Poly, Rational

z,t,u,k = symbols('z t u k')

def A_j(j):
    return sum(z**r/factorial(r) for r in range(1,j))

def C_j(j):
    return exp(z)-sum(z**r/factorial(r) for r in range(0,j))

def O_n(j,n):
    Aj,Cj=A_j(j),C_j(j)
    F=1/(1-Aj-t*Cj)
    return sp.expand(factorial(n)*sp.series(F,z,0,n+1).removeO().coeff(z,n))

def P_n(j,n):
    q=n//j
    On=O_n(j,n)
    return sp.expand((1-u)**q * On.subs(t,u/(1-u)))

def check_prompt_identity():
    print('== prompt identity ==')
    for j in range(2,7):
        Aj,Cj=A_j(j),C_j(j)
        diff=sp.simplify((1-Aj-t*Cj) - (1+t)*(1-Aj - t/(1+t)*exp(z)))
        print(f' j={j} diff={sp.simplify(diff)}')
        # correct
        diff2=sp.simplify((1-u)*(1-Aj-u/(1-u)*Cj) - (1-Aj*(1-2*u)-u*exp(z)))
        print(f'  correct={diff2}')

def Bnk(j,n):
    # B_{n,k} = n![z^n] e^{kz}/(1-A(1-2u))^(k+1) expanded in u? No, coupled.
    # Uncoupled wrong one: with (1-A-u e^z)
    Aj=A_j(j)
    res=0
    for m in range(0,n+1):
        dm = factorial(m)*sp.series((1-Aj)**(-(k+1)),z,0,m+1).removeO().coeff(z,m)
        dm = sp.expand(dm)
        res += binomial(n,m)*k**(n-m)*dm
    return sp.expand(res)

def check_wrong_decoupling():
    print('== wrong decoupling gives wrong P ==')
    for j in [3]:
        for n in range(0,9):
            q=n//j
            # wrong P from prompt identity
            # if identity were true: O_n = (1-u) sum u^k Bnk_wrong
            B=Bnk(j,n)
            # rational sum via falling basis
            Bpol=sp.expand(B)
            degB=sp.degree(Poly(Bpol,k))
            cps=[]
            for p in range(0,degB+1):
                cp=sum((-1)**(p-i)*binomial(p,i)*Bpol.subs(k,i) for i in range(0,p+1))
                cps.append(sp.simplify(cp))
            S=sum(cps[p]*u**p/(1-u)**(p+1) for p in range(0,len(cps)))
            Pwrong=sp.simplify((1-u)**(q+1)*S)
            Pright=P_n(j,n)
            print(f' j={j} n={n} equal={sp.simplify(Pwrong-Pright)==0}')

def check_correct_decoupling():
    print('== correct coupled decoupling ==')
    # F = (1-u)/(1-A(1-2u)-u e^z)
    # O_n = (1-u) n![z^n] sum_{k>=0} u^k e^{kz} (1-A(1-2u))^{-(k+1)}
    for j in [3,4]:
        Aj=A_j(j)
        for n in range(0,8):
            q=n//j
            # compute via series in z then u
            On=O_n(j,n)
            Pdirect=sp.expand((1-u)**q * On.subs(t,u/(1-u)))
            # coupled formula: coefficientwise in u
            # P = (1-u)^(q+1) * sum_{k>=0} u^k B_{n,k}(u) where B has u from (1-2u)
            # compute B_{n,k}(u) as polynomial in u,k
            # (1-A(1-2u))^{-(k+1)} = sum_p binom(k+p,p) A^p (1-2u)^p
            total=0
            for kk in range(0,n+5):
                Bk_u=0
                for m in range(0,n+1):
                    dm=0
                    Ap_pow=sp.Integer(1)
                    for p in range(0,m+1):
                        coeff=sp.Rational(factorial(m))*sp.series(Ap_pow,z,0,m+1).removeO().coeff(z,m)
                        coeff=sp.simplify(coeff)
                        if coeff!=0:
                            dm += coeff*binomial(kk+p,p)*(1-2*u)**p
                        Ap_pow=sp.expand(Ap_pow*Aj)
                    Bk_u += binomial(n,m)*kk**(n-m)*dm
                total += u**kk * sp.expand(Bk_u)
            Pcoupled = sp.expand((1-u)**(q+1)*total)
            # Pcoupled is a series in u; truncate to degree q
            Pu = Poly(Pcoupled,u)
            Ptrunc = sum(Pu.nth(i)*u**i for i in range(0,q+1))
            Ptrunc=sp.expand(Ptrunc)
            print(f' j={j} n={n} ok={sp.simplify(Ptrunc-Pdirect)==0}')

def check_realrootedness_u():
    print('== real-rootedness in u ==')
    for j in range(2,6):
        for n in range(0, min(20, 6*j)):
            P=P_n(j,n)
            q=n//j
            if q<=1:
                continue
            poly=sp.Poly(P,u)
            # count real roots by sturm
            from sympy import sturm, count_roots
            try:
                cr=count_roots(P,u)
            except Exception as e:
                cr=f'err:{e}'
            print(f' j={j} n={n} deg={q} real={cr}')

if __name__=='__main__':
    check_prompt_identity()
    check_wrong_decoupling()
    check_correct_decoupling()
    # check_realrootedness_u()  # can be slow
