exec(open('proof-final-threshold/family_probe.py').read().split('for j in [3,4,5]:')[0])
for j in [3,4,5]:
    N=16
    Oj=O(N,j); Hj=H(N,j)
    print('=== j',j)
    for n in range(j,N+1):
        fam=[('H',n,Hj[n])]+[('O',n-r,Oj[n-r]) for r in range(1,j) if n-r>=0]
        roots_fam=[roots(p) for _,_,p in fam]
        bad=[]
        for i in range(len(fam)):
            for k in range(i+1,len(fam)):
                if not interlace(roots_fam[i],roots_fam[k]):
                    bad.append((fam[i][:2],fam[k][:2]))
        if bad: print(' n',n,'bad',bad)
