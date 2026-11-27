exec(open('proof-final-threshold/family_probe.py').read().split('for j in [3,4,5]:')[0])
for j in [3,4,5]:
    N=18
    Oj=O(N,j); Hj=H(N,j)
    print('=== j',j)
    for m in range(0,N+1):
        bad=[]
        for tt in range(0,j):
            if m+tt>N: continue
            if not interlace(roots(Hj[m+tt]),roots(Oj[m])):
                bad.append(tt)
        if bad: print(' O_%d fails with H_{m+t} t='%m,bad)
