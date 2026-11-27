import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from scan import ordered_rows
from sympy import Poly, Symbol, Rational, solve, Matrix

t = Symbol('t')

def to_poly(coeffs):
    return sum(c*t**i for i,c in enumerate(coeffs))

for j in [2,3]:
    rows = [c for _, c in ordered_rows(j, 30)]
    # subsequence n = j*k + s for s in 0..j-1
    for s in range(j):
        seq = [rows[n] for n in range(s, len(rows), j) if n < len(rows)]
        # degrees
        degs = [len(c)-1 for c in seq]
        print('j',j,'s',s,'degs',degs[:8])
        # Try 3-term recurrence: seq[k+1] = (a_k*t + b_k)*seq[k] + c_k*seq[k-1]
        # solve for a_k,b_k,c_k using leading coefficients
        # Use k=2.. to have seq[k-1]
        oks=True
        for k in range(2, min(8, len(seq)-1)):
            pk = Poly(to_poly(seq[k]), t)
            pk1 = Poly(to_poly(seq[k+1]), t)
            pk_1 = Poly(to_poly(seq[k-1]), t)
            # leading coeff of pk1 = a_k * lc(pk) (if deg increases by 1)
            # assume deg increases by 1 each step
            lc1 = pk1.nth(pk1.degree())
            lc = pk.nth(pk.degree())
            a = Rational(lc1, lc)
            # constant term: pk1(0) = b_k*pk(0) + c_k*pk_1(0)
            # t coeff etc. Solve linear system using two evaluation points
            # Build equations by comparing coefficients
            # (a*t+b)*pk + c*pk_1 = pk1
            # Unknowns b,c (a known). Use coefficient of t^0 and t^1 maybe.
            import sympy as sp
            b,c = sp.symbols('b c')
            lhs = sp.expand((a*t + b)*to_poly(seq[k]) + c*to_poly(seq[k-1]))
            rhs = to_poly(seq[k+1])
            # pad
            dmax = max(sp.degree(lhs,t), sp.degree(rhs,t))
            eqs=[]
            for d in range(dmax+1):
                eqs.append(sp.Eq(sp.expand(lhs).coeff(t,d), sp.expand(rhs).coeff(t,d)))
            sol = sp.solve(eqs, [b,c], dict=True)
            if not sol:
                oks=False
                print('  no sol k',k)
                break
            else:
                # verify all coeffs
                sol=sol[0]
                lhs2 = sp.expand((a*t + sol[b])*to_poly(seq[k]) + sol[c]*to_poly(seq[k-1]))
                if lhs2 != rhs:
                    oks=False
                    print('  mismatch k',k, sol)
                    break
        print('  3-term recurrence holds:', oks)
