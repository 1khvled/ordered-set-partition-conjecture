import mpmath as mp
mp.mp.dps = 50
from saddle_explore import A_j, C_j, G_j, H_j, psi, alpha_pole

j=3
ap = alpha_pole(j)
print('ap', ap)
for a in [0.1,0.25,0.5,0.75,0.9]:
    mu = j*a/(1-a)
    lo = mp.mpf('1e-6')
    hi = ap - mp.mpf('1e-12')
    flo = psi(j,lo)-mu
    fhi = psi(j,hi)-mu
    print('a',a,'mu',mu,'flo',flo,'fhi',fhi)
    # also check middle
    for z in [0.1,0.3,0.5,0.7,0.73, ap-mp.mpf('1e-6'), ap-mp.mpf('1e-9')]:
        try:
            v = psi(j,mp.mpf(z))-mu
            print('  z',z,'f',v)
        except Exception as e:
            print('  z',z,'err',e)
