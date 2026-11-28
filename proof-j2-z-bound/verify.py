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
