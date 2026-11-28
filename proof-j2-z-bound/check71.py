from probe import probe, Qx, ev, dev
import numpy as np
for (n, av) in [(6,2),(8,2)]:
    R = Qx(n, av); V = Qx(n, av+1); P = Qx(n+2, av)
    for d in probe(n, av):
        s = d['s']
        # recompute Y
        from probe import Qx
        H = Qx(n-1, av+1)
        Y = 2*H + np.array([0.0]+list(dev(H)))
        lhs = ev(P,s)/av
        rhs = av*ev(R,s) + s*ev(dev(V),s) - n*s*ev(Y,s)
        print(n,av,"s=%.5f P/a=%.4f  RHS=%.4f  g=%.3f Y/R=%.3f"%(s,lhs,rhs,d['g'],d['YR']))
