from sympy import Rational, factorial, binomial, Poly, symbols, ZZ, QQ, div, resultant, discriminant, sqrt
from sympy import expand, cancel, series, Poly as P
import itertools

t,z = symbols('t z')

def A_j(j):
    return sum(z**r/factorial(r) for r in range(1,j))

def C_j(j):
    # C_j = e^z - sum_{r=0}^{j-1} z^r/r! as series; we only need polynomial coefficients up to N
    pass

# compute coefficients exactly up to N
def compute_H_O(j, N):
    # g_j = C_j/(1-A_j), H_j = 1/(1-t*g_j)
    # compute series coefficients in z up to N
    # Use series in QQ[t][[z]] represented as dict
    # power series as list of polynomials in t
    def series_add(a,b):
        L=max(len(a),len(b))
        res=[0]*L
        for i in range(L):
            res[i]=(a[i] if i<len(a) else 0) + (b[i] if i<len(b) else 0)
        return res
    def series_sub(a,b):
        L=max(len(a),len(b))
        res=[0]*L
        for i in range(L):
            res[i]=(a[i] if i<len(a) else 0) - (b[i] if i<len(b) else 0)
        return res
    def series_mul(a,b,N):
        res=[0]*(N+1)
        for i,ai in enumerate(a):
            if ai==0: continue
            for k,bk in enumerate(b):
                if bk==0: continue
                if i+k>N: break
                res[i+k]=res[i+k]+ai*bk
        return res
    def series_inv(a,N):
        # a[0] must be invertible constant
        assert a[0]!=0
        inv0 = 1/a[0]
        res=[0]*(N+1)
        res[0]=inv0
        for n in range(1,N+1):
            s=0
            for k in range(1,n+1):
                if k<len(a) and a[k]!=0 and res[n-k]!=0:
                    s += a[k]*res[n-k]
            res[n]= -inv0*s
        return res
    def series_scalar(a,c):
        return [c*x for x in a]

    # build A_j, C_j, D_j up to N
    A=[Rational(0)]*(N+1)
    for r in range(1,j):
        if r<=N: A[r]=Rational(1,factorial(r))
    one_minus_A=[Rational(0)]*(N+1)
    one_minus_A[0]=Rational(1)
    for r in range(1,j):
        if r<=N: one_minus_A[r]= -Rational(1,factorial(r))
    D=series_inv(one_minus_A,N)  # 1/(1-A_j)
    # C_j = e^z - E_{j-1}
    C=[Rational(0)]*(N+1)
    for r in range(j,N+1):
        C[r]=Rational(1,factorial(r))
    # g_j = C_j * D_j
    g=series_mul(C,D,N)
    # H_j(z,t) = 1/(1-t*g) as series in z with coeffs in Q[t]
    # compute iteratively: H = 1 + t*g*H
    H=[0]*(N+1)
    H[0]=1  # in Q[t]
    for n in range(1,N+1):
        # H[n] = t * sum_{k=1}^{n} g[k]*H[n-k]
        s=0
        for k in range(1,n+1):
            if g[k]!=0 and H[n-k]!=0:
                s = s + g[k]*H[n-k]
        H[n]=expand(t*s)
    # O_n = n! [z^n] D/(1-t*g) = n! [z^n] D*H
    DH=series_mul(D,H,N)
    O=[expand(factorial(n)*DH[n]) for n in range(N+1)]
    Hrows=[expand(factorial(n)*H[n]) for n in range(N+1)]  # H_{n,j}=n! h_{n,j}
    return g,D,Hrows,O

for j in [2,3,4]:
    N=12
    g,D,Hrows,O = compute_H_O(j,N)
    print(f'j={j}')
    for n in range(0,min(N+1,13)):
        hn=expand(Hrows[n])
        on=expand(O[n])
        degH = Poly(hn,t).degree() if hn!=0 else -1
        print(f' n={n}: H_n(0)={hn.subs(t,0)}, degH={degH}, O_n(0?) no, O_n has deg {Poly(on,t).degree()}')
        # check recurrence T1: O_n = H_n + sum_{r=1}^{j-1} C(n,r)*O_{n-r}
        if n>=0:
            rhs=hn
            for r in range(1,j):
                if n-r>=0:
                    rhs = rhs + binomial(n,r)*O[n-r]
            rhs=expand(rhs)
            if expand(on-rhs)!=0:
                print('  RECURRENCE FAIL', n)
    print('---')
