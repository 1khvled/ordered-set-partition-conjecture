import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-consecutive')
from verify_consecutive import ordered_rows, poly_from_coeffs
import sympy as sp

t=sp.Symbol('t')
O={}
for n,c in ordered_rows(2,10): O[n]=poly_from_coeffs(c)
# Compute W_n
def W(n):
    P=O[n]; Q=O[n+1]
    return sp.expand(P.as_expr()*sp.diff(Q.as_expr(),t)-sp.diff(P.as_expr(),t)*Q.as_expr())
for n in range(2,7):
    Wn=sp.Poly(W(n),t)
    Wn1=sp.Poly(W(n+1),t)
    # try divide W_{n+1} by W_n? 
    q,r=sp.div(Wn1.as_expr(), Wn.as_expr(), t)
    print(f"n={n} degWn={Wn.degree()} degWn1={Wn1.degree()} remainder_degree={sp.Poly(r,t).degree() if r!=0 else -1}")