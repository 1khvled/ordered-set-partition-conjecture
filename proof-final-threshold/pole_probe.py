import mpmath as mp
mp.mp.dps=40

def E(z,r):
    return sum(z**k/mp.factorial(k) for k in range(r+1))

def tau(z,j):
    N=2-E(z,j-1)
    C=mp.e**z-E(z,j-1)
    return N/C

for j in [2,3,4]:
    zs=[mp.mpf('0.7'),mp.mpf('1.3'),mp.mpf('2.1'),mp.mpf('3.4')]
    print('j',j)
    for z in zs:
        print(' ',float(z),'tau_{j-1}',mp.nstr(tau(z,j-1),12),'tau_j',mp.nstr(tau(z,j),12))
    # fit Moebius phi with phi(tau_{j-1}(z_i))=tau_j(z_i) using 3 pts, test 4th
    x=[tau(z,j-1) for z in zs]
    y=[tau(z,j) for z in zs]
    # solve a*x+b = y*(c*x+d), set d=1, unknown a,b,c
    A=mp.matrix([[x[k],1,-y[k]*x[k]] for k in range(3)])
    b=mp.matrix([y[k] for k in range(3)])
    sol=mp.lu_solve(A,b)
    a,bb,c=sol[0],sol[1],sol[2]
    pred=(a*x[3]+bb)/(c*x[3]+1) if abs(c*x[3]+1)>1e-30 else mp.inf
    print('  fitted a,b,c',mp.nstr(a,8),mp.nstr(bb,8),mp.nstr(c,8),'pred',mp.nstr(pred,12),'actual',mp.nstr(y[3] if False else tau(zs[3],j),12))
