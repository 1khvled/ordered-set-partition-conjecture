import sys, math
sys.path.insert(0, r"C:\Users\Abdelli\Documents\ordered-big-blocks\proof-final-inner")
from explore_H import inner_array
from fractions import Fraction

def H_coeffs(j,n,nmax=24):
    h=inner_array(j,nmax)
    row=h[n]; deg=n//j; fact=math.factorial(n)
    return [fact*row[m] for m in range(deg+1)]

for j in [2,3,4]:
    bad=0; total=0
    for n in range(1,24):
        cs=H_coeffs(j,n)
        # Newton inequalities for m=1..deg-1: a_m^2 >= (m+1)/m * (n? not degree) ... for polynomial degree d: c_m^2 >= (m+1)/m * (d-m+1)/(d-m) c_{m-1} c_{m+1}
        d=len(cs)-1
        for m in range(1,d):
            lhs=cs[m]*cs[m]
            rhs=Fraction(m+1,m)*Fraction(d-m+1,d-m)*cs[m-1]*cs[m+1]
            total+=1
            if lhs < rhs:
                bad+=1
                if bad<=5: print("FAIL",j,n,m,lhs,rhs)
    print(f"j={j} Newton bad {bad}/{total}")
