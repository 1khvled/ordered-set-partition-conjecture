import transform_lib as tl
import sympy as sp
t=sp.symbols('t')
j=3; N=14
O,H,D=tl.compute_O(j,N)
def roots(expr):
    if sp.degree(expr,t)<=0: return []
    return sorted([complex(x) for x in sp.nroots(expr)], key=lambda z:z.real)
def interlace(froots,groots):
    fr=sorted([z.real for z in froots]); gr=sorted([z.real for z in groots])
    if len(fr)==0 or len(gr)==0: return 'trivial'
    if abs(len(fr)-len(gr))>1: return False
    items=sorted([(z,'f') for z in fr]+[(z,'g') for z in gr])
    colors=[c for _,c in items]
    for i in range(1,len(colors)):
        if colors[i]==colors[i-1]: return False
    return True
bad=0; total=0
for a in range(0,N+1):
    for b in range(a+1,N+1):
        total+=1
        if not interlace(roots(O[a].as_expr()), roots(O[b].as_expr())):
            bad+=1
            if bad<=10: print('NONinterlace O',a,b)
print('j',j,'bad pairs among O',bad,'/',total)
