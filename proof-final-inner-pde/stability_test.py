import numpy as np
import cmath, math

def g(j,z):
    # C_j/N
    E=sum(z**r/math.factorial(r) for r in range(0,j))
    C=cmath.exp(z)-E
    A=sum(z**r/math.factorial(r) for r in range(1,j))
    N=1-A
    return C/N

for j in [2,3,4]:
    # sample upper half-plane
    bad=0; total=0
    for x in np.linspace(-3,3,41):
        for y in np.linspace(0.1,3,30):
            z=x+1j*y
            gz=g(j,z)
            total+=1
            # check imag of gz
            if gz.imag < -1e-9:
                bad+=1
    print(f"j={j}: samples with Im g<0: {bad}/{total}")
    # also check if 1/g in UHP for some z
    cnt=0
    for x in np.linspace(-3,3,41):
        for y in np.linspace(0.1,3,30):
            z=x+1j*y
            gz=g(j,z)
            if abs(gz)>1e-12:
                inv=1/gz
                if inv.imag>1e-9:
                    cnt+=1
    print(f"   z in UHP with Im(1/g)>0: {cnt}")
