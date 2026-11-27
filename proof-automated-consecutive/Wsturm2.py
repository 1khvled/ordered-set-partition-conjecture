import sys, math, time
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs
import sympy as sp

t = sp.Symbol('t')

def sturm_real_count(poly_expr):
    seq = sp.sturm(poly_expr, t)
    def eval_signs(val):
        nums=[]
        for s in seq:
            nums.append(sp.limit(s, t, val))
        signs=[0 if v==0 else (1 if v>0 else -1) for v in nums]
        nz=[s for s in signs if s!=0]
        ch=0
        for i in range(len(nz)-1):
            if nz[i]*nz[i+1]<0: ch+=1
        return ch
    return eval_signs(-sp.oo)-eval_signs(sp.oo)

def run(j,nmax):
    t0=time.time()
    O={}
    for n,coeffs in ordered_rows(j,nmax):
        O[n]=poly_from_coeffs(coeffs)
    bad=0; maxdeg=0
    for n in range(1,nmax):
        P=O[n]; Q=O[n+1]
        if P.degree()<=0: continue
        W=sp.expand(P.as_expr()*sp.diff(Q.as_expr(),t) - sp.diff(P.as_expr(),t)*Q.as_expr())
        maxdeg=max(maxdeg, sp.Poly(W,t).degree())
        nr=sturm_real_count(W)
        if nr!=0:
            bad+=1
            print(f"   BAD j={j} n={n} roots={nr}")
    print(f"j={j} nmax={nmax} maxdeg={maxdeg} bad={bad} time={time.time()-t0:.1f}s")

run(2,35)
run(3,25)