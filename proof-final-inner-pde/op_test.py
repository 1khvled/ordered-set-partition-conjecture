import numpy as np
import random

def roots_real(coeffs):
    # coeffs list low to high
    c = list(coeffs)
    while len(c)>1 and abs(c[-1])<1e-12: c.pop()
    if len(c)<=1: return True
    rts = np.roots(c[::-1])
    return max(abs(z.imag) for z in rts) < 1e-8

def random_realrooted(deg):
    # roots negative random
    rts = [-random.random()*5 - 0.1 for _ in range(deg)]
    # monic polynomial with those roots: prod (t - r)
    coeffs = [1.0]
    for r in rts:
        # multiply by (t - r): new = old shifted - r*old
        new = [0]*(len(coeffs)+1)
        for i,c in enumerate(coeffs):
            new[i+1]+=c
            new[i]-=r*c
        coeffs=new
    # low to high already? prod gives high-to-low? Let's ensure: coeffs[0] constant. Our recurrence produces constant (-r product) yes low-to-high.
    return coeffs

def T_op(p):
    # t * p'
    deg=len(p)-1
    out=[0]*(deg+1)  # degree deg? t*p' degree deg (if p degree deg, p' degree deg-1, times t degree deg)
    # p' coeff at k: (k+1) p[k+1]; times t -> index k+1
    for k in range(deg):
        out[k+1] = (k+1)*p[k+1]
    return out

def Q_op(p):
    # p + (t-1) p' = p - p' + t p'
    # compute p' low-high
    pp=[0]*(len(p)-1) if len(p)>1 else [0]
    for k in range(1,len(p)):
        pp[k-1]=k*p[k]
    # (t-1)p' = t p' - p'
    # t p': shift pp up by 1
    tp=[0]+pp
    # pad
    n=max(len(p),len(tp),len(pp))
    def pad(a): return a+[0]*(n-len(a))
    p=pad(p); tp=pad(tp); pp=pad(pp)
    return [p[i]-pp[i]+tp[i] for i in range(n)]

random.seed(0)
cnt=0; bad=0
for _ in range(2000):
    deg=random.randint(2,6)
    p=random_realrooted(deg)
    for op,name in [(T_op,'T'),(Q_op,'Q')]:
        q=op(p)
        if not roots_real(q):
            bad+=1
            if bad<=5:
                print("FAIL",name,p,q)
    cnt+=1
print("trials",cnt,"bad",bad)
