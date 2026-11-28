# `j=2` diagonal route after (H): adjacent step proved, cone still open

Date: 2026-11-28  
Scope: only `ordered-big-blocks/proof-j2-diagonal-after-h` was created/edited.  
No installs, builds, external agents, models, or browser use.

## 0. Outcome

**No all-`x` closure is obtained.** There is, however, a new rigorous positive result:

\[
\boxed{W\bigl(Q_{m,a},Q_{m-1,a+1}\bigr)(t)>0\quad(t\in\mathbb R)}
\tag{I1}
\]

for every nonconstant adjacent diagonal pair, where
\(Q_{m,a}=B_{m,a}/t\). Thus (H) is no longer the missing input for the
*adjacent* diagonal Wronskian. In fact (I1) follows from the stronger coupled
interlacing theorem in `proof-u-operator/NOTES.md` §7, the PF multiplier
\(\Lambda_a\), and the exact vertical identity; the coefficient inequality (H)
is not needed separately.

This does **not** close the diagonal-cone argument, for two exact reasons:

1. Equal-degree and terminal degree-zero components are not controlled by
   adjacent interlacing. The next nonautomatic diagonal inequality is
   \[
   \boxed{W\bigl(Q_{m,a},Q_{m-2,a+2}\bigr)(t)>0\quad(t\in\mathbb R).}
   \tag{I2}
   \]
   A rational feasible model satisfying (H), factorial-row real-rootedness
   and coupled interlacing, all component horizontal/vertical interlacings,
   and (I1), violates (I2).
2. The source component decomposition omits the \(m=0\) term \(n!x^n\).
   Consequently the final polynomial is
   \[
   O_n(x,t)=n!x^n+tS_{n,x}(t),
   \]
   not merely a positive combination of the \(Q_{m,a}\). Even a compatible
   nonconstant diagonal family does not control the terminal constants or the
   top term. An exact feasible model satisfies (I1) and (I2), yet its \(S\)
   already fails real-rootedness at \(x=1\); a simpler independent example has
   real-rooted \(S\) but non-real-rooted \(1+tS\).

Accordingly, the exact remaining route conditions are (I2) for the
nonconstant diagonal blocks **and** the boundary inequality (J) in §6. The
all-`x` \(j=2\) conjecture remains open. No consecutive-row nonvanishing or
interlacing statement (R) is obtained either.

## 1. Definitions and a necessary boundary correction

Let

\[
A_{m,r}={m\brace r}_{\ge2},\qquad
G_m(t)=\frac{E_m(t)}t=\sum_{r\ge1}r!A_{m,r}t^{r-1},
\]

and

\[
Q_{m,a}(t)=\frac{B_{m,a}(t)}t
 =\sum_{r\ge1}(a)_rA_{m,r}t^{r-1}.
\]

For fixed \(M\), write

\[
D_a^{(M)}(t)=Q_{M-a,a}(t),\qquad 1\le a\le M-2.
\]

The source `proof-j2-resultants/NOTES.md` §3.1 formula
\(P_n^{(k)}=(n!/m!)B_{m,a}\), \(m=n-k,\ a=k+1\), is correct for \(m\ge1\)
(with \(B_{1,a}=0\)), but not for \(m=0\): the \(r=0\) term in
\((1-tC)^{-a}\) contributes \(n!\) when \(k=n\). The exact decomposition is

\[
\boxed{
O_n(x,t)=n!x^n+
\sum_{k=0}^{n-2}x^k\frac{n!}{(n-k)!}B_{n-k,k+1}(t).
}
\tag{1}
\]

Equivalently, with \(M=n+1\),

\[
O_n(x,t)=n!x^n+tS_{n,x}(t),\qquad
S_{n,x}(t)=\sum_{a=1}^{M-2}
\frac{n!}{(M-a)!}x^{a-1}D_a^{(M)}(t).
\tag{2}
\]

`audit.py` checks (1) exactly for \(2\le n\le12\); every residual is zero.
This top term is not optional for the all-`x` implication.

## 2. What (H) gives on the diagonal

Put \(m=M-a\) and

\[
c_{a,r}=[t^{r-1}]D_a^{(M)}=(a)_rA_{m,r},\qquad
R_r=\frac{c_{a+1,r}}{c_{a,r}}.
\]

For \(r\ge2\), inequality (H), in the form

\[
(r-1)A_{m,r}A_{m-1,r-1}
 \ge rA_{m,r-1}A_{m-1,r},
\]

gives

\[
\frac{R_r}{R_{r-1}}
 =\frac{a+r}{a+r-1}
  \frac{A_{m-1,r}A_{m,r-1}}
       {A_{m-1,r-1}A_{m,r}}
 \le \frac{(a+r)(r-1)}{r(a+r-1)}<1.
\tag{3}
\]

For \(r=1\), (H) with \(k=1\) gives
\(A_{m,2}\ge2A_{m-1,2}\), hence

\[
\frac{R_1}{R_0}
 =\frac{a+2}{a+1}\frac{A_{m-1,2}}{A_{m,2}}
 \le\frac{a+2}{2(a+1)}<1.
\tag{4}
\]

Thus (H) is exactly strict diagonal TP2 in coefficient-ratio form:

\[
R_0>R_1>R_2>\cdots.
\tag{5}
\]

**Direction sample.** For \(M=8,a=1\), the actual coefficients are

\[
D_1=1+112t+630t^2,\qquad D_2=2+150t+360t^2,
\]

so \((R_0,R_1,R_2)=(2,75/56,4/7)\), strictly decreasing. At \(r=2\),

\[
R_2/R_1=32/75<3/4,
\]

which is the upper bound in (3). The direction is therefore
“cross-ratio \(\le\), factor \(<1\), hence \(R_r<R_{r-1}\).”

Coefficient TP2 alone is not a cone proof; §5 gives exact feasible failures.

## 3. The adjacent diagonal step is proved

### 3.1 The PF multiplier

For a polynomial \(f(t)=\sum_{j\ge0}f_jt^j\), define

\[
(\Lambda_af)(t)=\sum_{j\ge0}\binom{a+j}{j}f_jt^j.
\]

Since

\[
\sum_{j\ge0}\binom{a+j}{j}z^j=(1-z)^{-a},
\]

\(\Lambda_a\) is a Pólya-frequency multiplier sequence. It preserves
real-rootedness and strict interlacing; all multipliers are positive, so no
active coefficient or degree is lost. Directly from the coefficients,

\[
\boxed{Q_{m,a}=\Lambda_aG_m.}
\tag{6}
\]

Both \(\Lambda_a\) and \(\mathcal U=2+t\,d/dt\) are diagonal in the monomial
basis, so

\[
\boxed{\mathcal UQ_{m,a}=\Lambda_a\mathcal UG_m.}
\tag{7}
\]

The coefficient residual in (7) is
\((j+2)\binom{a+j}{j}f_j-(j+2)\binom{a+j}{j}f_j=0\).

### 3.2 Apply the coupled factorial-row theorem

For \(m\ge5\), `proof-u-operator/NOTES.md` §7 proves, in the same orientation,

\[
G_{m-1}\preceq G_m,\qquad
\mathcal UG_{m-1}\preceq G_m.
\tag{8}
\]

Applying \(\Lambda_a\) to both pairs and using (7) gives

\[
Q_{m-1,a}\preceq Q_{m,a},\qquad
\mathcal UQ_{m-1,a}\preceq Q_{m,a}.
\tag{9}
\]

With the Wronskian convention \(W(p,q)=p'q-pq'\), the orientation in (9) is
\(W(p,q)>0\) on \(\mathbb R\).

### 3.3 Use the vertical identity

Dividing Theorem B by \(t\) gives

\[
Q_{m,a+1}
 =Q_{m,a}+\frac1a(tQ_{m,a})'
 =\frac1a\bigl(\mathcal UQ_{m,a}+(a-1)Q_{m,a}\bigr).
\tag{10}
\]

Set \(p=Q_{m,a}\), \(q=Q_{m-1,a}\). By (9), both \(q\) and \(\mathcal Uq\)
interlace \(p\) in the same orientation. For \(a>1\), (10) is a positive
combination of these two compatible polynomials; for \(a=1\), it is simply
\(\mathcal Uq\). The two-polynomial Obreschkoff cone therefore gives

\[
Q_{m-1,a+1}\preceq Q_{m,a},
\]

which is exactly (I1).

**Direction sample.** At \(M=7,a=1\),

\[
p=1+50t+90t^2,\qquad q=2+60t,
\]

and

\[
W(p,q)=40+360t+5400t^2>0,
\]

because its discriminant is \(360^2-4(5400)(40)=-734400<0\) and its leading
coefficient is positive.

### 3.4 Small and inactive cases

* \(m=2,3\): \(Q_{m,a}\) is constant, so its Wronskian with another constant
  is zero; these pairs are omitted from strict (I1).
* \(m=4\): \(Q_{4,a}\) is linear and \(Q_{3,a+1}\) is constant, so
  \(W(Q_{4,a},Q_{3,a+1})=Q_{4,a}'Q_{3,a+1}>0\).
* \(m=1\): \(B_{1,a}=0\).
* \(m=0\): the separate top term \(n!x^n\) in (1).
* \(a=1\): the \((a-1)q\) term in (10) is zero and is handled directly.
* Inactive \(A_{m,r}=0\) entries extend all coefficient inequalities
  trivially; they are not used as positive active coefficients.

## 4. Why adjacent interlacing is not the full diagonal cone

The degree pattern is

\[
\deg D_a^{(M)}=\left\lfloor\frac{M-a}{2}\right\rfloor-1.
\]

Consecutive degrees differ by \(0\) or \(1\), and every equal-degree block has
at most two members. However:

* terminal \(D\)-polynomials are constants;
* a constant added to a higher-degree polynomial is not an interlacing
  partner in the Obreschkoff sense;
* equal-degree blocks can be adjacent to a higher-degree block on one side
  and constants on the other.

Therefore (I1) is not sufficient for all positive combinations in (2).
The first additional nonautomatic Wronskian condition for a nonconstant
triple is (I2). Since
\(\deg D_a-\deg D_{a+2}=1\), its required direction is unambiguously

\[
D_a^{(M)}\ \text{surrounds}\ D_{a+2}^{(M)}
\quad\Longleftrightarrow\quad
W(D_a^{(M)},D_{a+2}^{(M)})>0.
\tag{11}
\]

Even (I1)+(I2) does not control the terminal constants or the top term; §6
shows this by an exact feasible model.

## 5. Feasible falsification: (H)+factorial rows+(I1) do not imply (I2)

The following rational rows are a **logical countermodel**, not a claim that
the associated-Stirling recurrence holds:

\[
\begin{array}{c|ccc}
m&A_{m,1}&A_{m,2}&A_{m,3}\\ \hline
4&1&1/100&0\\
5&1&3/10&0\\
6&1&3/4&1/20
\end{array}
\tag{12}
\]

All active coefficients are positive. The exact (H) residuals are

\[
\begin{array}{c|c}
(m,k)&kA_{m+1,k+1}A_{m,k}-(k+1)A_{m+1,k}A_{m,k+1}\\ \hline
(4,1)&7/25\\
(4,2)&0\\
(5,1)&3/20\\
(5,2)&3/100
\end{array}
\]

so every residual is nonnegative, with the required direction.

The factorial polynomials are

\[
G_4=1+\frac{t}{50},\qquad
G_5=1+\frac{3t}{5},\qquad
G_6=1+\frac{3t}{2}+\frac{3t^2}{10}.
\]

The roots of \(G_6\) are \((-15\pm\sqrt{105})/6\). The pairs required by the
coupled factorial-row theorem satisfy, exactly,

\[
W(G_5,G_4)>0,\quad W(G_5,\mathcal UG_4)>0,\quad
W(G_6,G_5)>0,\quad W(G_6,\mathcal UG_5)>0
\]

on \(\mathbb R\). The unweighted Bóna–Mező horizontal pairs and the
\(Q_{6,a}\)-versus-\(Q_{5,a}\) pairs for \(a=1,2,3\) also have strictly positive
Wronskians on \(\mathbb R\). The vertical identity (10) has zero residual for
every tested \(m,a\).

For \(M=7\), the first three diagonal elements are

\[
D_1=1+\frac32t+\frac3{10}t^2,\qquad
D_2=2+\frac95t,\qquad
D_3=3+\frac3{25}t.
\]

Both adjacent inequalities hold:

\[
W(D_1,D_2)=\frac65+\frac65t+\frac{27}{50}t^2>0,
\qquad
W(D_2,D_3)=\frac{129}{25}>0.
\]

But the distance-two inequality fails:

\[
W(D_1,D_3)=\frac{219}{50}+\frac95t+\frac9{250}t^2
\]

has two negative real zeros (its discriminant is \(8154/3125>0\)), so it is
not positive on \(\mathbb R\). Consequently,

\[
\operatorname{Disc}_t(D_1+D_3)=-\frac{5439}{2500}<0.
\tag{13}
\]

This is a feasible violation of the proposed implication using only the
listed coefficient/root constraints. It does not contradict the actual fixed
associated-Stirling triangle. It proves that a successful argument must use
additional structure, essentially the associated-Stirling recurrence or a
direct proof of (I2) for the actual \(A_{m,r}\).

## 6. The terminal/top boundary inequality is independent

### 6.1 Exact boundary condition

Let \(d=\lfloor n/2\rfloor\). Suppose first that \(S=S_{n,x}\) is real-rooted
with negative roots. Then \(H(t)=tS(t)\) has the roots of \(S\) together with
\(0\). For \(d\ge2\), define

\[
\mathcal C(S)=
\min\bigl\{-H(\tau):H'(\tau)=0,\ \tau<0,\ H(\tau)<0\bigr\},
\qquad H(t)=tS(t).
\]

The minimum is over the negative local minima of \(H\); values tending to
zero at roots of \(H\) are not admissible limiting critical levels. The exact
condition for \(n!x^n+tS(t)\) to have \(d\) simple negative roots is

\[
\boxed{n!x^n<\mathcal C(S_{n,x}).}
\tag{J}
\]

For \(d=1\), the polynomial is linear and (J) is automatic. Equality in (J)
is a double-root boundary; larger \(n!x^n\) shifts the graph upward and loses
two real roots.

**Direction audit.** For \(S=(t+1)(t+2)\), one has

\[
(tS)'=2+6t+3t^2.
\]

On \((-1,0)\), the relevant critical point is
\(t_0=-1+1/\sqrt3\), and

\[
(tS)''(t_0)=6+6t_0>0,\qquad
t_0S(t_0)=-\frac{2}{3\sqrt3}.
\]

Thus \(\mathcal C(S)=2/(3\sqrt3)\). Taking the top coefficient \(c=1\)
violates \(c<\mathcal C(S)\), and indeed

\[
1+tS=t^3+3t^2+2t+1
\]

has discriminant \(-23<0\). This checks the monotonicity direction: raising
\(c\) raises the graph and decreases the real-root count.

### 6.2 A model satisfying (I1)+(I2) still fails the boundary

Change only \(A_{4,2}\) in (12) from \(1/100\) to \(1/10\). The (H) residuals
become \(1/10,0,3/20,3/100\), and all factorial, unweighted, horizontal, and
vertical checks still pass. Now

\[
D_1=1+\frac32t+\frac3{10}t^2,\qquad
D_2=2+\frac95t,\qquad
D_3=3+\frac65t,
\]

and

\[
W(D_1,D_3)=\frac{33}{10}+\frac95t+\frac9{25}t^2>0
\]

on \(\mathbb R\) (discriminant \(-189/125<0\)). Hence both (I1) and (I2) hold.

For \(n=6\), \(D_4=4\), \(D_5=5\), and (2), after dividing by \(6!\), is

\[
\frac{O_6(x,t)}{6!}
 =\frac{x^6}{t}+D_1+6xD_2+30x^2D_3+480x^3+1800x^4.
\]

At \(x=1\),

\[
S=2383+\frac{483}{10}t+\frac3{10}t^2,
\]

but

\[
\operatorname{Disc}_t S
 =\left(\frac{483}{10}\right)^2
  -4(2383)\frac3{10}
 =-\frac{52671}{100}<0.
\tag{14}
\]

Thus the terminal constants already destroy real-rootedness of \(S\), before
the top \(x^6/t\) term is added. Adding it also fails:

\[
10(1+tS)=3t^3+483t^2+23830t+10,
\]

whose discriminant is

\[
-29908514659080<0.
\tag{15}
\]

This is an exact feasible counterexample to the assertion that (H), the
factorial-row results, component interlacing, (I1), and even (I2) suffice for
the all-`x` conclusion.

As with (12), this model does not satisfy the associated-Stirling recurrence.
It refutes the stated logical implication, not the actual fixed triangle.

## 7. Source-example audit

The nontransitivity example in `proof-j2-resultants/NOTES.md` §4.3 is not a
valid strict-interlacing example as written:

\[
f=(t+1)(t+4),\qquad g=(t+2)(t+3)
\]

have merged root labels \(f,g,g,f\), not \(f,g,f,g\) or \(g,f,g,f\).
Thus that particular pair does not strictly interlace. The conclusion that
interlacing is not transitive remains true, but the displayed witness must be
replaced. The countermodels in §§5–6 supply exact replacements for the
specific implications needed here.

## 8. Policy audit and finite evidence

* **Direction:** all Wronskians required in (I1), (I2), and (11) are \(>0\),
  not \(<0\); explicit samples and counterexamples are given above.
* **Monotonicity:** (3)–(5) print the exact discrete ratio and its \(<1\)
  direction. For (J), the derivative, critical point, second derivative, and
  endpoint/critical-value comparison are printed in §6.1.
* **Feasible falsification:** (12) violates (I2) while satisfying every listed
  input; the \(A_{4,2}=1/10\) variant violates the boundary while satisfying
  (I1)+(I2). These are exact rational witnesses, not numerical searches.
* **Boundary/degeneracy:** \(m=0,1,2,3\), \(a=1\), inactive indices, equal
  degrees, terminal constants, and the top \(x^n\) term are covered in §§1,
  3.4, and 6.
* **Symbolic residuals:** `audit.py` prints zero for the corrected component
  decomposition, the vertical identity, and the inactive (H) cases.
* **Finite checks are not proof:** on the actual triangle, exact Sturm checks
  find no failure of (H) through \(m=17\), no failure of (I1) through \(M=18\),
  and no failure of (I2) through \(M=18\). This is finite evidence only.

## 9. Exact remaining work

For the actual associated-Stirling triangle, the unresolved inequalities are

\[
\boxed{W(Q_{m,a},Q_{m-2,a+2})>0\quad\text{on }\mathbb R}
\tag{I2}
\]

for all nonautomatic \(m,a\), together with the terminal/top boundary
condition (J), or a stronger direct coefficient/Bezoutian inequality that
implies both. Neither is proved by (H) plus factorial-row real-rootedness
alone. Therefore the all-`x` \(j=2\) conjecture is not closed by this route.

## 10. Reproduction

From this directory:

```powershell
py -3.14 audit.py
```

The script uses only the Python standard library and exact rational arithmetic.
