exec(open('proof-final-threshold/family_probe.py').read().split('for j in [3,4,5]:')[0])
for j in [2,3,4,5]:
    N=20
    Hj=H(N,j)
    print('=== j',j)
    bad=[]
    for a in range(j,N+1):
        if not Hj[a].equals(0):
            for b in range(a+1,N+1):
                if Hj[b].equals(0): continue
                if not interlace(roots(Hj[a]),roots(Hj[b])):
                    bad.append((a,b))
    print('bad pairs',bad[:20], 'count',len(bad))
