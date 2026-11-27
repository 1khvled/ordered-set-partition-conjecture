# Corrected WRS reattack for the actual ordered big-block family

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-wrs-corrected` was created or edited.

## Outcome

**No unconditional all-index proof of WRS was obtained.**  The alleged
`j=2, n=63` counterexample is retracted, and the retraction is reproduced
below by exact rational arithmetic.

The reattack does produce three useful exact results:

1. WRS is equivalent to positivity of a Wronskian on the whole negative axis,
   not merely to positivity of its coefficients.
2. The full first-block history gives the exact identity
   \[
   W(O_n,O_{n-j})
   =\sum_{r=1}^{j-1}\binom nr W(O_{n-r},O_{n-j})
    +W(H_{n,j},O_{n-j}),
   \]
   where \(H_{n,j}\) is the inner chapter row.  Thus the entire older tail is
   one concrete cross-family Wronskian.
3. Two tempting strengthenings are blocked by exact actual-family or formal
   calculations: the natural all-lag weighted inequality fails at an actual
   consecutive pair, and coefficientwise Wronskian positivity is insufficient
   even for a formal recurrent example.

The precise remaining inequality is displayed in Section 8.  It is a
quantitative cross-family residue/Bezoutian bound; neither ordinary inner
interlacing nor the positive lag-`j` coefficient-minor recurrence supplies it.

---

## 1. Family, recurrence, and the correct target

Fix \(j\ge2\).  Write

\[
B_j(z,t)=\sum_{r=1}^{j-1}\frac{z^r}{r!}
+t\left(e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}\right),
\qquad
F_j(z,t)=\frac1{1-B_j(z,t)},
\]

\[
O_n=n![z^n]F_j(z,t),\qquad O_0=1.
\]

The actual first-block recurrence is

\[
O_n=\sum_{r=1}^{j-1}\binom nr O_{n-r}
+t\sum_{r=j}^{n}\binom nr O_{n-r}.
\tag{1}
\]

Put \(T=t\partial_t\), \(Q(p)=p+(t-1)p'\), and

\[
A_n=nO_{n-1}+2TO_{n-1}
+\sum_{r=2}^{j-1}\binom nr O_{n-r},
\qquad c=\binom{n-1}{j-1}.
\]

The differential recurrence is

\[
O_n=A_n+c\,tQ(O_{n-j}).
\tag{2}
\]

For \(q=O_{n-j}\) and a root \(\mu\) of \(q\),

\[
tQ(q)(\mu)=\mu(\mu-1)q'(\mu).
\]

Hence WRS is exactly

\[
-\frac{A_n(\mu)}{q'(\mu)}>c\mu(\mu-1)
\quad\Longleftrightarrow\quad
\frac{O_n(\mu)}{q'(\mu)}<0.
\tag{3}
\]

No factor of \(\mu\) belongs in the denominator \(q'(\mu)\).

---

## 2. Wronskian form of WRS

Use the orientation

\[
W(p,q)=p'q-pq'.
\tag{4}
\]

Assume the standing induction hypothesis that \(q\) has degree \(e\), simple
negative roots \(\mu_1,\ldots,\mu_e\), and \(p=O_n\) has degree \(e+1\) with
positive leading coefficient.  The partial-fraction expansion is

\[
\frac{p(t)}{q(t)}
=at+b+\sum_{i=1}^{e}\frac{r_i}{t-\mu_i},
\qquad
r_i=\frac{p(\mu_i)}{q'(\mu_i)},\quad a>0.
\]

Therefore

\[
W(p,q)=q^2\left(\frac pq\right)'
=q^2\left(a-\sum_i\frac{r_i}{(t-\mu_i)^2}\right).
\tag{5}
\]

If every \(r_i<0\), the parenthesis in (5) is strictly positive for every real
\(t\).  Conversely, if \(W(p,q)>0\) for every real \(t\), then near each pole
\(\mu_i\),

\[
W(p,q)\sim-\frac{r_iq'(\mu_i)^2}{(t-\mu_i)^2},
\]

so \(r_i<0\).  Thus, under the simple-real-root hypothesis,

\[
\boxed{\text{WRS}\quad\Longleftrightarrow\quad W(O_n,O_{n-j})(t)>0
\text{ for every }t\in\mathbb R.}
\tag{6}
\]

This equivalence explains why coefficientwise positivity of \(W\) is not
enough: it proves positivity for \(t\ge0\), while WRS needs the negative axis.

---

## 3. Exact retraction audit of `j=2, n=63`

The old script formed the derivative coefficient list as

```text
[0, q_1, 2q_2, 3q_3, ...]
```

and then evaluated that list as an ascending-coefficient polynomial.  This is
\(xq'(x)\), not \(q'(x)\).  The correct list is

```text
[q_1, 2q_2, 3q_3, ...].
```

`audit_corrected.py` reproduces the exact certificate at

\[
x_0=-0.301049864539336319020771115593462774725514449468323282248633,
\qquad \varepsilon=10^{-55}.
\]

It proves by an exact rational Sturm sequence that \(q=O_{61,2}\) has exactly
one root in \((x_0-\varepsilon,x_0+\varepsilon)\).  Exact Horner evaluation
with a derivative remainder bound then gives

\[
q'(x_0)>0,\qquad O_{63,2}(x_0)<0.
\]

Consequently

\[
\frac{O_{63,2}(\mu)}{q'(\mu)}<0,
\]

so WRS holds at the certified root.  The residue sign is \(-1\), not \(+1\).
For comparison only, a 200-digit evaluation with the corrected derivative
gives the older-tail margin

\[
M(\mu)=4.972349660011042767890309503057084255795\ldots>0.
\]

That decimal is not used as a certificate.  **There is no actual-family WRS
counterexample at `j=2,n=63`.**

---

## 4. Full actual history and the exact Wronskian tail

Define the inner row

\[
H_{n,j}(t)=n![z^n]\frac1{1-tg_j(z)},
\qquad
g_j(z)=\frac{C_j(z)}{1-A_j(z)}.
\]

The exact structural identities are

\[
H_{n,j}=t\sum_{r=j}^{n}\binom nr O_{n-r},
\tag{7}
\]

and therefore

\[
O_n=H_{n,j}+\sum_{r=1}^{j-1}\binom nr O_{n-r}.
\tag{8}
\]

Let \(q=O_{n-j}\).  Applying the bilinear map \(p\mapsto W(p,q)\) to (8)
gives the full-history identity

\[
\boxed{
W(O_n,q)=
\sum_{r=1}^{j-1}\binom nr W(O_{n-r},q)
+W(H_{n,j},q).
}
\tag{9}
\]

Thus the whole signed older tail is not an unspecified collection of terms:
it is exactly the cross-family Wronskian \(W(H_{n,j},O_{n-j})\).

A second exact identity is essential when expanding the tail directly:

\[
\boxed{W(tp,q)=pq+tW(p,q).}
\tag{10}
\]

The product \(pq\) in (10) cannot be dropped.  For example, with actual rows
\(q=O_{4,2}=24+45t+6t^2\) and \(p=O_{0,2}=1\),

\[
W(t,q)=q-tq'=24-6t^2,
\]

which has a negative coefficient.  This term occurs in the first-block
expansion of \(O_{6,2}\).  Hence a proof that treats every older summand as
just \(tW(O_m,q)\) is algebraically wrong, and a termwise positivity proof
fails even before evaluating at negative \(t\).

For the same actual pair \(n=6,j=2\),

\[
H_{6,2}=t(517+830t+90t^2),
\]

\[
W(H_{6,2},O_{4,2})
=12408+39840t+40728t^2+8100t^3+540t^4,
\]

and

\[
W(O_{6,2},O_{4,2})
=24792+62880t+59232t^2+8100t^3+540t^4.
\]

The recent term \(6W(O_{5,2},O_{4,2})\) supplies exactly the difference.
This illustrates the cancellation that a valid proof must control.

---

## 5. What the positive minor recurrence proves

The lag-`j` adjacent coefficient minor identity is

\[
\begin{aligned}
\Delta^{(j)}_{n,k}={}&
(n+2k)\Delta^{(j-1)}_{n-1,k}
+2u_{k+1}q_k\\
&+\sum_{r=2}^{j-1}\binom nr\Delta^{(j-r)}_{n-r,k}\\
&+c\bigl((k+1)q_k^2-q_kq_{k+1}-kq_{k-1}q_{k+1}\bigr),
\end{aligned}
\tag{11}
\]

where \(q=O_{n-j}\) and \(u=O_{n-1}\).  The actual first-block recurrence gives
\(u_{k+1}\ge cq_{k+1}\), and log-concavity of \(q\) makes the last line strictly
positive on the active coefficient range.  This is a genuine positive
recurrence for those coefficient minors.

It is not a residue bound.  In Wronskian language, coefficient minors control
formal coefficient ratios near \(t=0\); WRS controls the signs of residues at
negative roots, equivalently positivity on the whole negative axis by (6).
The formal recurrent obstruction in `proof-final-wrs/NOTES.md` remains valid
as a logical obstruction, but it is **not** an actual-family counterexample.

A short coefficientwise illustration is enough.  For the same formal data

\[
q=t+1,\qquad
r=17t^2+\frac{973}{25}t+\frac{552}{25},
\]

one has

\[
W(r,q)=17t^2+34t+\frac{421}{25},
\]

whose coefficients are all positive, but

\[
W(r,q)(-1)=-\frac4{25}<0.
\]

Thus even coefficientwise Wronskian positivity does not imply WRS.

---

## 6. A natural all-lag strengthening fails on an actual pair

For \(1\le d\le j\), consider the stronger weighted statement

\[
-\frac{O_n(\mu)}{O_{n-d}'(\mu)}
>\binom{n-1}{d-1}\mu(\mu-1),
\qquad O_{n-d}(\mu)=0.
\tag{12}
\]

The desired WRS is the case \(d=j\).  However, (12) already fails for an
actual consecutive pair, \(j=2,n=4,d=1\).  Here

\[
q=O_{3,2}=6+7t,\qquad \mu=-\frac67,
\qquad p=O_{4,2}=24+45t+6t^2.
\]

Exact evaluation gives

\[
-\frac{p(\mu)}{q'(\mu)}=\frac{498}{343},
\qquad
\mu(\mu-1)=\frac{78}{49}=\frac{546}{343}.
\]

Thus (12) fails by \(48/343\).  This is **not** a WRS counterexample: the
target for \(j=2\) uses \(d=j=2\), not \(d=1\).  It only shows that the
natural binomial-weighted all-lag invariant cannot be used as the missing
induction hypothesis.

---

## 7. Matrix/transform audit

Two simple total-positivity shortcuts also fail.

### Outer convolution matrix

The outer binomial convolution is represented, after positive diagonal
scaling, by the ordinary Riordan array

\[
K^{(j)}=(D_j(z),z),\qquad D_j(z)=\frac1{1-A_j(z)}.
\]

For \(j=3\), the ordinary coefficients of \(D_3\) begin

\[
1,\quad 1,\quad \frac32,\quad 2,\quad\ldots
\]

and therefore

\[
\det K^{(3)}[\{1,2\},\{0,1\}]
=
\det\begin{pmatrix}1&1\\[2pt]3/2&1\end{pmatrix}
=-\frac12.
\tag{13}
\]

So the outer matrix is not TN for \(j\ge3\).  A proof of the form “both
factors are TN, hence their product is TN” is unavailable.

### Natural inner production matrix

For the inner coefficient array \(h_{n,m}=[z^n]g_j(z)^m\), the natural
production matrix is the Toeplitz matrix of the ordinary coefficients of
\(g_j=C_jD_j\).  For \(j=3\),

\[
g_3=\frac16,\qquad g_4=\frac5{24},\qquad g_5=\frac3{10},
\]

and its minor with rows \(\{4,5\}\), columns \(\{0,1\}\) is

\[
g_4^2-g_3g_5=\frac{25}{576}-\frac1{20}
=-\frac{19}{2880}.
\tag{14}
\]

Thus the natural inner production matrix is not TN.  This does not assert
that the inner array itself is non-TN; it blocks that particular production
matrix proof.  The published associated-Bell/interlacing theorem controls
inner rows, but it does not bound the cross-family term
\(W(H_{n,j},O_{n-j})\) in (9).

---

## 8. Precise remaining gap

Let

\[
q=O_{n-j},\qquad
\mathcal W_{n,j}(t)=
\sum_{r=1}^{j-1}\binom nr W(O_{n-r},q)(t)
+W(H_{n,j},q)(t).
\]

By (9), \(\mathcal W_{n,j}=W(O_n,q)\).  The exact all-index statement still
needed is

\[
\boxed{
\mathcal W_{n,j}(t)>0
\quad\text{for every }j\ge2,\ n\ge j,\text{ and }t\le0.
}
\tag{15}
\]

By (6), (15) is equivalent to WRS under the simple-real-root induction
hypothesis.

A weaker coefficientwise target is to prove

\[
[t^k]W(H_{n,j},O_{n-j})\ge0
\quad\text{for all }k,
\tag{16}
\]

and simultaneously propagate coefficientwise positivity of the recent
Wronskians in (9).  Exact checks support both statements, but even if they are
proved they yield positivity only for \(t\ge0\).  The missing content is the
quantitative negative-axis/Bezoutian domination in (15), using the actual
chapter history rather than only local compatibility.

No derivation of (15), or of the equivalent main-part residue lower bound

\[
-\frac{A_n(\mu)}{O_{n-j}'(\mu)}
>\binom{n-1}{j-1}\mu(\mu-1),
\]

was obtained in this track.

---

## 9. Reproducible exact checks

Run from this directory:

```powershell
py -3.14 audit_corrected.py
py -3.14 check_invariants.py
```

`audit_corrected.py` uses FLINT only to obtain rigorous isolating Arb balls.
Every reported root sign is then rechecked by exact rational Horner evaluation
with a derivative remainder bound.  It certifies:

- the retracted `j=2,n=63` root has `q'>0`, `O_63<0`, hence negative residue;
- WRS root signs for `j=2..6` through `n=40,36,32,30,30`, respectively:
  361, 176, 91, 55, and 40 root checks, all passing.

`check_invariants.py` additionally verifies with integer arithmetic:

- the product-rule identity (10) and the explicit negative-coefficient term;
- the full-history identity (9) on the displayed `j=2,n=6` example;
- coefficientwise Wronskian signs for every row pair through
  `(j,nmax)=(2,60),(3,45),(4,36),(5,30),(6,30)`;
- coefficientwise nonnegativity of `W(H_{n,j},O_{n-j})` over the same ranges;
- the exact actual-pair obstruction to the all-lag strengthening (12);
- the exact matrix minors (13) and (14).

These finite checks are evidence and regression certificates, not an
all-index proof.  No packages were installed, no builds were run, and no
large files or caches were created.
