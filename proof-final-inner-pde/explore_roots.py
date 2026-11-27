import sys, math
sys.path.insert(0, r"C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-inner")
from explore_H import inner_array
from fractions import Fraction

def poly_coeffs(j,n):
    h = inner_array(j, 20)
    # h[n][m] already h/n!?? inner_array stores h[n][m] = [t^m] H_n? Let's inspect: doc says h[n][m] = [t^m] H_n? Actually returns h[n][m] as Fraction; explore prints. We'll just use as coefficient of t^m in H_n (not divided? Let's check).
    row = h[n]
    # degree
    deg = n//j
    coeffs = [float(row[m]) for m in range(deg+1)]
    return coeffs

if __name__=="__main__":
    for j in [2,3,4]:
        print("j=",j)
        for n in range(0,16):
            cs = poly_coeffs(j,n)
            # polynomial sum cs[m] t^m. Use numpy roots
            import numpy as np
            # remove trailing zeros
            while len(cs)>1 and abs(cs[-1])<1e-15: cs.pop()
            if len(cs)<=1:
                print(f"  n={n} coeffs={cs}")
                continue
            rts = np.roots(cs[::-1])  # highest degree first
            imag = [abs(z.imag) for z in rts]
            print(f"  n={n} deg={len(cs)-1} max|imag|={max(imag):.2e} roots={sorted(rts,key=lambda z:(z.real,z.imag))}")
