"""Exact symbolic residuals and direction checks for proof-j2-z-bound."""
import sympy as sp
from ops import a, x, A, Q

def sep(s):
    print("\n" + "="*70); print(s); print("="*70)

def alpha(n, av):
    al=[]; v=1
    for j in range(n//2+1):
        if j>0: v*=av+j
        al.append(sp.Integer(v)*A[n][j+1] if j+1<len(A[n]) else sp.Integer(0))
    return al

sep("V1: coefficient ratios r_j strictly increasing (Newton, exact samples)")
bad=0; tot=0
for n in range(6,25):
    for av in [1,2,3,5,10]:
        al=alpha(n,av)
        for j in range(len(al)-2):
            if al[j+2]:
                tot+=1
                bad+=int(not (al[j+1]**2 > al[j]*al[j+2]))
print(f"r_j increasing: checked={tot} violations={bad}")

sep("V2: mu=(a+1)A_{n,2}^2/((a+2)A_{n,3})>6 exactly, n<=30,a<=100")
bad=0; tot=0
for n in range(6,31):
    for av in range(1,101):
        mu=sp.Rational(av+1)*A[n][2]**2/(sp.Rational(av+2)*A[n][3])
        tot+=1
        bad+=int(not (mu>6))
print(f"mu>6: checked={tot} violations={bad}")
print("min mu at n=6,a=1:", sp.Rational(2)*A[6][2]**2/(sp.Rational(3)*A[6][3]))

sep("V3: recurrences A_{m+1,2}=2A_{m,2}+m ; A_{m+1,3}=3A_{m,3}+mA_{m-1,2}")
bad=0
for m in range(3,30):
    bad+=int(A[m+1][2]!=2*A[m][2]+m)
    bad+=int(A[m+1][3]!=3*A[m][3]+m*A[m-1][2])
print("recurrences bad:", bad)

sep("V4: p(t)<0 at t=r1/4 (sample exact)")
for (n,av) in [(6,2),(8,2),(12,5)]:
    al=alpha(n,av)[:n//2]
    r1=al[1]/al[2]; mu=al[1]**2/al[2]; t=r1/4
    pval=sum((-1)**j*al[j]*t**j for j in range(len(al)))
    print(f"  n={n} a={av}: p(t)={float(pval):.6g} <= 1-mu/6={float(1-mu/6):.6g} <0")

sep("V5: exact (M),(V) residuals")
bad=0
for n in range(4,12):
    for av in [1,2,3]:
        q=Q(n,av,x); qm=Q(n-1,av,x); qv=Q(n-2,av+1,x); qa1=Q(n,av+1,x)
        bad+=int(sp.expand(q-(sp.diff(x*qm,x)+(n-1)*av*x*qv))!=0)
        bad+=int(sp.expand(av*qa1-((av+1)*q+x*sp.diff(q,x)))!=0)
print("residual bad:", bad)

sep("V6: direction 2a+1 < 3a+2 for a>0")
for av in [1,2,3,5,10]:
    print(f"  a={av}: (a+1)^2-a^2={2*av+1} < 3a+2={3*av+2}: {2*av+1<3*av+2}")

sep("V7: exact rational obstruction to a generic non-first-root extension")
av = 20
rr = sp.expand((1 - 2*x) * (1 - sp.Rational(13, 14)*x)
               * (1 - sp.Rational(10, 11)*x))
vv = sp.expand((av + 1)*rr + x*sp.diff(rr, x)) / av
rhos = [sp.Rational(1, 2), sp.Rational(14, 13), sp.Rational(11, 10)]
us = [sp.simplify(1/(1 - rho)) for rho in rhos]
z1 = sp.simplify(sp.diff(rr, x, 2).subs(x, 1) / rr.subs(x, 1))
print("R =", sp.factor(rr))
print("V =", sp.factor(vv))
print("V(1) =", sp.simplify(vv.subs(x, 1)))
print("u_j =", us, "sum =", sum(us))
print("Z(1) =", z1, "threshold =", 3*av + 2,
      "violation:", bool(z1 >= 3*av + 2))
print("signs at 0,rho1,1,rho2,rho3 =",
      [sp.sign(vv.subs(x, z)) for z in [0, rhos[0], 1, rhos[1], rhos[2]]])
print("degree(V) =", sp.Poly(vv, x).degree())

sep("V8: n=4,a=2 direct smallest-root direction sample")
nn, aa = 4, 2
rr = Q(nn, aa, -x)
vv = Q(nn, aa + 1, -x)
hh = Q(nn - 1, aa + 1, -x)
kk = Q(nn - 2, aa + 2, -x)
yy = sp.expand(2*hh + x*sp.diff(hh, x))
sigma = sp.Rational(1, 12)
rho = sp.Rational(1, 9)
p1 = sp.simplify(sigma/(rho - sigma))
z0 = sp.simplify(sigma**2*sp.diff(rr, x, 2).subs(x, sigma)
                 / rr.subs(x, sigma))
print("R =", rr, "V =", vv, "H =", hh, "K =", kk, "Y =", yy)
print("sigma =", sigma, "rho =", rho, "V(sigma) =", vv.subs(x, sigma))
print("p1 =", p1, "> a:", bool(p1 > aa), "Z =", z0)
print("M-root residual =", sp.simplify(
    (hh + sigma*sp.diff(hh, x)).subs(x, sigma)
    - (nn - 1)*(aa + 1)*sigma*kk.subs(x, sigma)))
print("Y/2 =", sp.simplify(yy.subs(x, sigma)/2))
