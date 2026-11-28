# C1 cross term for 2-associated Stirling rows: retraction and exact gap

Date: 2026-11-28
Scope: only `ordered-big-blocks/proof-associated-stirling-c1` was created/edited.

## Status

**The all-index C1 claim is retracted.** The previous proof had a monotonicity sign
error. The two prior-row SLC inequalities do imply the two constraints derived below,
but those constraints alone do not imply $C1\ge0$. No repaired all-index proof was
obtained in this directory.

This note records the valid identities, the corrected optimization, an exact formal
obstruction, and the remaining gap.

## 1. Notation and valid consequences of prior SLC

Let $A_{r,k}=\{r\brace k\}_{\ge2}$ and
$$
A_{r,k}=kA_{r-1,k}+(r-1)A_{r-2,k-1}.
$$
Put $n=m-1$ and
$$
A=A_{n,k},\quad B=A_{n,k-1},\quad C=A_{n,k+1},\quad
P=A_{n-1,k-1},\quad Q=A_{n-1,k-2},\quad R=A_{n-1,k}.
$$
The target is
$$
C1(m,k)=2k^3AP-(k+1)^2\bigl((k-1)BR+(k+1)CQ\bigr).
$$
Assume, as in the proposed SLC induction,
$$
SLC(n,k)=k^2A^2-(k+1)^2BC\ge0, \tag{1}
$$
$$
SLC(n-1,k-1)=(k-1)^2P^2-k^2QR\ge0. \tag{2}
$$
For $k\ge2$, $A,B,P,R>0$. Set
$$
X=AP,\qquad Y=BR,\qquad Z=CQ,\qquad
y=\frac YX,\qquad z=\frac ZX,
$$
and
$$
c=\frac{k-1}{k+1},\qquad r=\frac{k-1}{k},\qquad
\alpha=(k+1)^2(k-1),\qquad \beta=(k+1)^3.
$$

### 1.1 Product constraint

The exact identity
$$
(k+1)^2BCQR-(k-1)^2A^2P^2
=-QR\,SLC(n,k)-A^2\,SLC(n-1,k-1)
$$
gives
$$
yz\le c^2. \tag{3}
$$

### 1.2 Consecutive-row ratio constraint

The recurrence gives
$$
A=kR+(n-1)P,\qquad B=(k-1)P+(n-1)Q,
$$
hence
$$
(k-1)AP-kBR=(n-1)\bigl((k-1)P^2-kQR\bigr).
$$
Also
$$
(k-1)P^2-kQR
=\frac{SLC(n-1,k-1)}{k-1}+\frac{k}{k-1}QR\ge0.
$$
Therefore
$$
0<y\le r. \tag{4}
$$
Equivalently,
$$
\frac{A_{m-1,k}A_{m-2,k-1}}
{A_{m-1,k-1}A_{m-2,k}}\ge\frac{k}{k-1}.
$$
This is a valid lower bound on the cross-row ratio $X/Y$, not the lower bound on
$y=Y/X$ that the failed optimization needed.

## 2. Corrected optimization and exact obstruction

The normalized target is
$$
\frac{C1}{X}=2k^3-\alpha y-\beta z.
$$
From (3), $0\le z\le c^2/y$, so
$$
\frac{C1}{X}\ge 2k^3-f(y),\qquad
f(y)=\alpha y+\frac{\beta c^2}{y}. \tag{5}
$$
However,
$$
f'(y)=\alpha-\frac{\beta c^2}{y^2}<0
\quad\text{on }0<y\le r,
$$
because
$$
\beta c^2-\alpha r^2=\frac{(k+1)(k-1)^2}{k^2}>0.
$$
Thus the correct direction is
$$
f(y)\ge f(r), \tag{6}
$$
not $f(y)\le f(r)$. The endpoint is
$$
f(r)=2k^3-3k+\frac1k.
$$
Consequently the old computation only shows
$$
2k^3-f(y)\le 3k-\frac1k,
$$
which is an upper bound on the lower estimate in (5), not a lower bound on $C1/X$.

The constraints (3)–(4) alone are insufficient. For any $k\ge2$, choose arbitrarily
small $y>0$ and set $z=c^2/y$. Then (3)–(4) hold, while
$$
\frac{C1}{X}=2k^3-\alpha y-\frac{\beta c^2}{y}\longrightarrow-\infty
\qquad(y\to0^+). \tag{7}
$$
This is a formal obstruction to proving C1 from only those two constraints; it is not
a counterexample using actual associated Stirling numbers.

For the requested $k=2$ check,
$$
y=\frac1{10},\qquad z=\frac1{30}
$$
satisfies
$$
yz=\frac1{300}\le\frac19=c^2,\qquad y\le\frac12=r,
$$
but
$$
\frac{C1}{X}=16-9\left(\frac1{10}\right)-27\left(\frac1{30}\right)
=\frac{71}{5}=14.2.
$$
Thus the retracted value $11/2=5.5$ is neither a valid lower nor a valid upper bound
on $C1/X$.

A negative feasible constraint example is obtained at $k=2$ by taking
$y=1/100$ and $z=100/9$. Then $yz=1/9$ and
$$
\frac{C1}{X}=16-\frac9{100}-27\frac{100}{9}
=-\frac{28409}{100}<0.
$$
Again, this is a formal point satisfying only (3)–(4), not an actual Stirling-row
configuration.

## 3. Valid partial cases

### 3.1 $k=1$

Here $B=P=Q=0$, so $C1(m,1)=0$.

### 3.2 $k=2$

Here $Q=A_{n-1,0}=0$, so $Z=0$. Using (4),
$$
\frac{C1}{X}=2k^3-\alpha y
\ge2k^3-\alpha r
=k^3+2k-\frac1k>0.
$$
At $k=2$ this gives $C1/X\ge23/2$.

In fact the $k=2$ case is explicit. For $s\ge2$,
$$
A_{s,2}=2^{s-1}-s-1,
$$
because partitions into two blocks are counted by $(2^s-2)/2$, and exactly $s$ of
them have a singleton block. Therefore, with $n=m-1\ge5$,
$$
C1(m,2)=16A_{n,2}-9A_{n-1,2}
=23\cdot2^{n-2}-7n-16>0.
$$

### 3.3 Top boundary $m=2k+2$

Here $n=2k+1$ and $C=A_{n,k+1}=0$, hence $Z=0$. Again (4) gives
$$
\frac{C1}{X}\ge k^3+2k-\frac1k>0.
$$
The boundary SLC inequality (1) is trivially true because $C=0$.

## 4. Exact remaining gap

The unresolved range is
$$
k\ge3,\qquad m\ge2k+3,
$$
where $Z=CQ>0$. To repair the proof by the present optimization, one needs a
rigorous additional constraint on $(y,z)$, such as a lower bound on $y$ strong enough
that
$$
\alpha y+\frac{\beta c^2}{y}\le2k^3
$$
throughout the actual feasible range, or a sharper upper bound on $z$ than $c^2/y$.
For precision, the product constraint alone would suffice if one could prove
$$
y\ge L_-(k):=
\frac{k^3-\sqrt{3k^4-3k^2+1}}{(k+1)^2(k-1)},
$$
because $L_-(k)$ is the smaller root of
$\alpha y^2-2k^3y+\beta c^2=0$. No such lower bound, or substitute joint constraint,
has been proved here.

Since C1 remains unproved on that range, the SLC induction is **not** closed by this
note. The earlier claim that C0, C1, and C2 were all controlled is withdrawn.

## 5. Symbolic verification

`verify_c1.py` checks the valid identities, the corrected direction (6), the formal
limit (7), and the requested $k=2$ constraint example. It performs no finite
enumeration and makes no all-index C1 claim.

```powershell
py -3.14 verify_c1.py
```
