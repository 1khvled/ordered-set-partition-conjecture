import sympy as sp
t = sp.symbols('t')

def rows(j, nmax):
    O = [sp.Integer(1)]
    for n in range(1, nmax+1):
        s = sp.Integer(0)
        for r in range(1, min(j-1, n)+1):
            s += sp.binomial(n, r) * O[n-r]
        for r in range(j, n+1):
            s += t * sp.binomial(n, r) * O[n-r]
        O.append(sp.expand(s))
    return O

def W(p,q):
    return sp.expand(sp.diff(p,t)*q - p*sp.diff(q,t))

def H_n(O,j,n):
    s = sp.Integer(0)
    for r in range(j,n+1):
        s += sp.binomial(n,r)*O[n-r]
    return sp.expand(t*s)

for j in range(2,8):
    nmax = min(5*j, 30)
    O = rows(j,nmax)
    bad = 0
    tested = 0
    for n in range(j, nmax+1):
        for d in range(1, j+1):
            h = H_n(O,j,n)
            w = W(h, O[n-d])
            if w == 0: 
                continue
            tested += 1
            rc = sp.count_roots(sp.Poly(w,t))
            if rc != 0:
                bad += 1
                print('BAD', j,n,d,rc)
    print(f'j={j} tested={tested} bad={bad}')
