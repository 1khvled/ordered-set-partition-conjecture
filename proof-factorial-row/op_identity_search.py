"""Search for identity linking W(U p, V R) to W(p,R) and monotone quantities."""
import sympy as sp
from explore import T
xs = sp.symbols('x')
for r in range(6, 10):
    p = T[r-2].as_expr(); q = T[r-3].as_expr(); R = xs*q
    Up = 2*p + xs*sp.diff(p,xs)
    VR = sp.diff(xs*R, xs)
    Wnew = sp.expand(sp.diff(Up,xs)*VR - Up*sp.diff(VR,xs))
    Wold = sp.expand(sp.diff(p,xs)*R - p*sp.diff(R,xs))
    # candidate: Wnew = c * Wold * (something polynomial)
    ratio = sp.factor(sp.cancel(Wnew/Wold))
    print('r=%d ratio Wnew/Wold =' % r, ratio)
