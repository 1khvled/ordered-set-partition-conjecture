"""Search exact identity for P/R at roots of V in terms of Z and coupled quantities."""
import sympy as sp
from ops import a, x, R, Q

for n in [6,8]:
    for av in [2,3]:
        Rn=R(n,av); Vn=R(n,av+1); Hn=R(n-1,av+1); Sn=R(n+1,av); Pn=R(n+2,av)
        Yn=sp.expand(2*Hn+x*sp.diff(Hn,x))
        # P = (xS)' - (n+1)a x V  ; S = (xR)' - n a x H
        S2=sp.expand(sp.diff(x*Rn,x)-n*av*x*Hn)
        P2=sp.expand(sp.diff(x*S2,x)-(n+1)*av*x*Vn)
        print(f"n={n} a={av}: P identity residual =", sp.simplify(Pn-P2))
        # At root of V: claim P = x^2 R'' - (3a+2)R - n a (x H + Y - H)
        claim=sp.expand(x**2*sp.diff(Rn,x,2)-(3*av+2)*Rn-n*av*(x*Hn+Yn-Hn))
        rem=sp.rem(sp.expand(Pn-claim),Vn,x)
        print("   rem(P-claim mod V) =", sp.factor(rem))
