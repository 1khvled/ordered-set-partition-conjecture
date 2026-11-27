import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs
import sympy as sp

t = sp.Symbol('t')

def sturm_real_roots(poly_expr):
    p=sp.Poly(poly_expr,t,domain='ZZ')
    # number of distinct real roots via Sturm
    sqf = sp.sqf_list(poly_expr)
    # use sturm
    seq = sp.sturm(poly_expr, t)
    # count sign changes at -inf and +inf
    def signs_at(val):
        vals=[sp.limit(s, t, val) if val in (sp.oo,-sp.oo) else s.subs(t,val) for s in seq]
        # evaluate limits
        nums=[]
        for s in seq:
            if val==sp.oo:
                nums.append(sp.limit(s,t,sp.oo))
            else:
                nums.append(sp.limit(s,t,-sp.oo))
        # sign
        signs=[]
        for v in nums:
            if v==0:
                signs.append(0)
            else:
                signs.append(1 if v>0 else -1)
        # remove zeros
        nz=[s for s in signs if s!=0]
        changes=0
        for i in range(len(nz)-1):
            if nz[i]*nz[i+1]<0:
                changes+=1
        return changes
    cpos=signs_at(sp.oo)
    cneg=signs_at(-sp.oo)
    return cneg-cpos

def Wsturm(j,nmax):
    O={}
    for n,coeffs in ordered_rows(j,nmax):
        O[n]=poly_from_coeffs(coeffs)
    bad=0
    for n in range(1,nmax):
        P=O[n]; Q=O[n+1]
        if P.degree()<=0: continue
        W=sp.expand(P.as_expr()*sp.diff(Q.as_expr(),t) - sp.diff(P.as_expr(),t)*Q.as_expr())
        nr=sturm_real_roots(W)
        if nr!=0:
            bad+=1
            print(f"   BAD j={j} n={n} realroots={nr}")
    print(f"j={j} nmax={nmax} bad={bad}")

for j in [2,3,4,5]:
    Wsturm(j,25 if j==2 else 20)