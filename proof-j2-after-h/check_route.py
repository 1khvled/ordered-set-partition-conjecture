from fractions import Fraction
import math

# 2-associated Stirling numbers of 2nd kind: A_{m,k} = {m brace k}_{>=2}
# recurrence: A_{m,k} = k*A_{m-1,k} + (m-1)*A_{m-2,k-1}
# base: A_{0,0}=1, A_{m,0}=0 for m>0, A_{m,k}=0 if m<2k or k<0

def compute_A(max_m):
    A = [[Fraction(0) for _ in range(max_m+1)] for _ in range(max_m+1)]
    A[0][0] = Fraction(1)
    for m in range(1, max_m+1):
        for k in range(0, max_m+1):
            if m < 2*k:
                A[m][k] = Fraction(0)
            elif k == 0:
                A[m][k] = Fraction(0)
            else:
                val = Fraction(0)
                if m-1 >= 0:
                    val += k * A[m-1][k]
                if m-2 >= 0 and k-1 >= 0:
                    val += (m-1) * A[m-2][k-1]
                A[m][k] = val
    return A

# O_n(t) coefficients: [t^k]O_n = n! * sum_{m=2k}^{n} (n-m+1)_k / m! * A_{m,k}

def falling(a, k):
    r = Fraction(1)
    for j in range(k):
        r *= (a-j)
    return r

def O_coeffs(n, A):
    coeffs = [Fraction(0)] * (n//2 + 1)
    nfac = math.factorial(n)
    for k in range(0, n//2 + 1):
        s = Fraction(0)
        for m in range(2*k, n+1):
            w = Fraction(nfac, math.factorial(m)) * falling(n-m+1, k)
            s += w * A[m][k]
        coeffs[k] = s
    return coeffs

# Stirling numbers of 2nd kind S2(k,i)

def S2(max_k):
    S = [[Fraction(0) for _ in range(max_k+1)] for _ in range(max_k+1)]
    S[0][0] = Fraction(1)
    for k in range(1, max_k+1):
        for i in range(1, k+1):
            S[k][i] = i*S[k-1][i] + S[k-1][i-1]
    return S

max_n = 14
A = compute_A(max_n)
S = S2(max_n)

# check ratio monotonicity c_{n+1,i}/c_{n,i} increasing in i
print('ratio monotonicity c_{n+1,i}/c_{n,i} increasing in i:')
for n in range(2, max_n-1):
    d1 = O_coeffs(n, A)
    d2 = O_coeffs(n+1, A)
    c1 = [Fraction(0)]*(len(d1))
    c2 = [Fraction(0)]*(len(d2))
    for i in range(len(d1)):
        s=Fraction(0)
        for k in range(i,len(d1)):
            s += S[k][i]*d1[k]
        c1[i]=s
    for i in range(len(d2)):
        s=Fraction(0)
        for k in range(i,len(d2)):
            s += S[k][i]*d2[k]
        c2[i]=s
    ok = True
    for i in range(len(c1)-1):
        lhs = c2[i+1]*c1[i] - c2[i]*c1[i+1]
        if lhs < 0:
            ok = False
            print(f'  FAIL n={n} i={i} lhs={lhs}')
    print(f'  n={n}: {ok}')

# check TP2 of d_{n,k}=[t^k]O_n
print('\nTP2 of d_{n,k}=[t^k]O_n:')
for n in range(2, max_n-1):
    d1 = O_coeffs(n, A)
    d2 = O_coeffs(n+1, A)
    m = min(len(d1), len(d2)-1)
    ok = True
    for k in range(m-1):
        lhs = d1[k]*d2[k+1] - d1[k+1]*d2[k]
        if lhs < 0:
            ok = False
            print(f'  FAIL n={n} k={k} lhs={lhs}')
    print(f'  n={n}: {ok}')

# check TP2 of c_{n,i}
print('\nTP2 of c_{n,i} (falling factorial coeffs):')
for n in range(2, max_n-1):
    d1 = O_coeffs(n, A)
    d2 = O_coeffs(n+1, A)
    c1 = [Fraction(0)]*(len(d1))
    c2 = [Fraction(0)]*(len(d2))
    for i in range(len(d1)):
        s=Fraction(0)
        for k in range(i,len(d1)): s += S[k][i]*d1[k]
        c1[i]=s
    for i in range(len(d2)):
        s=Fraction(0)
        for k in range(i,len(d2)): s += S[k][i]*d2[k]
        c2[i]=s
    m = min(len(c1), len(c2)-1)
    ok = True
    for i in range(m-1):
        lhs = c1[i]*c2[i+1] - c1[i+1]*c2[i]
        if lhs < 0:
            ok = False
            print(f'  FAIL n={n} i={i} lhs={lhs}')
    print(f'  n={n}: {ok}')
