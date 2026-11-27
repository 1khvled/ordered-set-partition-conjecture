import mpmath as mp
mp.mp.dps = 30
from saddle_explore import A_j, C_j, G_j, H_j, psi, alpha_pole

for j in [3]:
    ap = alpha_pole(j)
    print('j', j, 'ap', ap)
    for z in [0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.73]:
        z = mp.mpf(str(z))
        print(' z', z, 'psi', psi(j, z), 'H', H_j(j, z))
