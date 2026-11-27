exec(open('proof-final-threshold/family_probe.py').read().split('for j in [3,4,5]:')[0])
for j in [3,4,5]:
    N=18
    Oj=O(N,j); Hj=H(N,j)
    print('=== j',j)
    worst=0
    for m in range(0,N+1):
        for tt in range(0,N-m+1):
            if Hj[m+tt].equals(0): continue
            if not interlace(roots(Hj[m+tt]),roots(Oj[m])):
                if tt>worst: worst=tt
                print(' first fail O_%d with t=%d'%(m,tt)); break
    print(' max failing t',worst)
