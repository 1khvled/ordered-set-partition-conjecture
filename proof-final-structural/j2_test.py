import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from scan import ordered_rows
from fractions import Fraction

# Simpler family H_r(t) = r! [z^r] 1/(1 - t*g_2(z)), g_2 = (e^z-1-z)/(1-z)
# Compute via recurrence: 1/(1-tg) = sum h_r z^r/r! ? Let's just compute power series coeffs (ordinary) Hr_coeffs
# and then h_r polynomial = r! * coeff (as poly in t).

# We'll compute coefficients of 1/(1 - t*g) as power series in z with polynomial coeffs in t (sympy).
import sympy as sp
z, t = sp.symbols('z t')
# g_2 series up to degree 12
N=12
# g coeffs from earlier
def g_coeffs(j,N):
    # reuse
    import importlib.util
    spec=importlib.util.spec_from_file_location('pc', r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-structural\pf_check.py')
    # just recompute inline
    pass

# compute g series via rational
from fractions import Fraction
def series_exp(N):
    return [Fraction(1, math.factorial(r)) for r in range(N+1)]
def series_inv_one_minus_A(j,N):
    A=[Fraction(0)]*(N+1)
    for r in range(1,j):
        if r<=N: A[r]=Fraction(1, math.factorial(r))
    D=[Fraction(0)]*(N+1); D[0]=Fraction(1)
    for n in range(1,N+1):
        s=Fraction(0)
        for r in range(1,min(j,n+1)):
            if r<=j-1: s+=A[r]*D[n-r]
        D[n]=s
    return D
def mul_series(a,b,N):
    c=[Fraction(0)]*(N+1)
    for i in range(N+1):
        for k in range(N+1-i): c[i+k]+=a[i]*b[k]
    return c
j=2
N=10
C=[Fraction(0)]*(N+1)
exp=series_exp(N)
for r in range(N+1):
    if r>=j: C[r]=exp[r]
D=series_inv_one_minus_A(j,N)
g=mul_series(C,D,N)

# Now H = 1/(1 - t*g) = sum_{m>=0} t^m g^m. Coeff [z^r] = polynomial in t.
# Compute by series inversion: H*(1 - t*g)=1. Let H coeffs Hr (polynomials in t) as lists.
# H_0 = 1. H_n = t * sum_{k=1}^{n} g_k * H_{n-k} (since (1-tg)H=1 => H = 1 + t g H)
# g_k are Fractions.
H=[None]*(N+1)
H[0]=sp.Poly(1, t)
for n in range(1,N+1):
    s=sp.Poly(0,t)
    for k in range(1,n+1):
        if k < len(g):
            s = s + sp.Poly(t*g[k]*H[n-k].as_expr(), t)  # g[k] Fraction -> ok
    # H_n = s (as poly) but need polynomial in t with rational coeffs
    H[n]=sp.Poly(sp.simplify(s.as_expr()), t)
    # multiply? H_n should be polynomial in t degree floor(n/j). Let's just keep.

# Check real-rootedness of h_r = r! * H_r? Actually H_r = [z^r]1/(1-tg) (ordinary coeff). 
# The "row polynomial" in Theorem 4.4 is r_n(t) = [x^n] 1/(1-t(f-f0)) (ordinary). So r_n = H_n (as poly in t). Check roots.
for r in range(0, N+1):
    poly = sp.Poly(H[r].as_expr(), t)
    roots = [complex(root) for root in sp.nroots(poly.as_expr())] if poly.degree()>0 else []
    nonreal = [rt for rt in roots if abs(rt.imag)>1e-9]
    print('r',r,'deg',poly.degree(),'nonreal_count',len(nonreal))
    if nonreal:
        print('   nonreal', nonreal[:3])
