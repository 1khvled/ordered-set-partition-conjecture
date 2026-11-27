import transform_lib as tl
import sympy as sp
t=sp.symbols('t')
j=3; N=10
O,H,D=tl.compute_O(j,N)
def roots(expr):
    if sp.degree(expr,t)<=0: return []
    return sorted([complex(x) for x in sp.nroots(expr)], key=lambda z:z.real)
def interlace(froots,groots):
    # check if real roots alternate: merge sorted; require no two same? 
    fr=sorted([z.real for z in froots]); gr=sorted([z.real for z in groots])
    if len(fr)==0 or len(gr)==0: return 'trivial'
    # degree diff <=1
    if abs(len(fr)-len(gr))>1: return False
    # interlacing: for sorted combined, colors alternate
    items=sorted([(z,'f') for z in fr]+[(z,'g') for z in gr])
    colors=[c for _,c in items]
    # allow equal? no
    for i in range(1,len(colors)):
        if colors[i]==colors[i-1]: 
            return False
    return True
for n in [7,8,9,10]:
    hf=roots(H[n].as_expr())
    for r in range(1,j):
        of=roots(O[n-r].as_expr())
        print('n',n,'r',r,'H-O interlace?', interlace(hf,of), 'degH',len(hf),'degO',len(of))
