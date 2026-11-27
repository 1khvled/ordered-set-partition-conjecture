# Weighted root-separation notes for the original ordered big-block family

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-weighted` was created or edited.

## 1. Family and conventions

Fix `j >= 2` and use the **original** series

\[
 B_j(z,t)=\sum_{r=1}^{j-1}\frac{z^r}{r!}
 +t\left(e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}\right),
 \qquad F_j(z,t)=\frac1{1-B_j(z,t)}.
\]

Thus `B_j(0,t)=0`, `F_j(0,t)=1`, and

\[
 O_n=n![z^n]F_j(z,t),\qquad O_0=1.
\]

No normalization correction is used anywhere in this note.

Put

\[
 E_{j-1}(z)=\sum_{r=0}^{j-1}\frac{z^r}{r!},\quad
 C_j(z)=e^z-E_{j-1}(z),\quad
 w_j(z)=\frac{z^{j-1}}{(j-1)!},\quad T=t\partial_t.
\]

Then `C_j'=C_j+w_j`. If `A_j=E_{j-1}-1`, then
`A_j'=A_j+1-w_j`, and hence

\[
 B_j'=B_j+1+(t-1)w_j.
\]

Since `F_j'=B_j'F_j^2` and `B_jF_j^2=F_j^2-F_j`,

\[
 \boxed{F_j'=\bigl(2+(t-1)w_j\bigr)F_j^2-F_j.} \tag{1}
\]

Also `(F_j)_t=C_jF_j^2`. These are the original-family Riccati identities.

The universal recurrence is

\[
 \boxed{
 O_n=R_n(O_{n-1})+S_n+C_{n,j}\,tQ(O_{n-j}),} \tag{2}
\]

where

\[
 R_n(p)=np+2Tp,\qquad
 S_n=\sum_{r=2}^{j-1}\binom nr O_{n-r},\qquad
 C_{n,j}=\binom{n-1}{j-1},\qquad
 Q(p)=p+(t-1)p'.
\]

Negative-index rows and the last term for `n<j` are omitted. The sum is empty for `j=2`.

For `p(t)=sum_k p_k t^k`, the signed coefficient operator is exactly

\[
 \boxed{tQ(p)=\sum_{k\ge1}k(p_{k-1}-p_k)t^k.} \tag{3}
\]

This follows directly from `tQ(p)=tp+t(t-1)p'`. It is the coefficient-level source of the obstruction.

---

## 2. The compatible-window sign lemma

The following standard lemma is useful here.

**Lemma 1.** Let a finite family of real-rooted polynomials with positive leading coefficients lie in one common compatible window, with the usual orientation: older equal-degree rows have their roots to the left, and a row of degree one less has one root in each gap of the newer row. If `rho` is a root of the newest member `p`, then every older member `q` satisfies

\[
 \operatorname{sgn}q(\rho)=\operatorname{sgn}p'(\rho). \tag{4}
\]

When `p` is newer than `q` and `rho` is a root of `q`, one also has

\[
 \operatorname{sgn}p(\rho)=-\operatorname{sgn}q'(\rho), \tag{5}
\]

both when `deg p=deg q+1` and when the degrees are equal with the newer row shifted right.

**Proof.** In a common window, the numbers of roots of `p` and `q` to the right of `rho` have the same parity in (4). Multiplying the signs contributed by those roots proves (4). For (5), if `q` has roots `mu_i` and `rho=mu_i`, then `q'(rho)` has sign `(-1)^(e-i)`, where `e=deg q`. If `p` has one higher degree, it has one additional root to the right. If the degrees are equal and `p` is shifted right, `rho` lies one gap earlier relative to the roots of `p`. In either case the parity is opposite. `Q.E.D.`

---

## 3. Important correction to the proposed IVT argument

At a root `mu` of `q=O_{n-j}`, the last term in (2) generally does **not** vanish. Instead,

\[
 \boxed{
 tQ(q)(\mu)=\mu(\mu-1)q'(\mu),} \tag{6}
\]

which is nonzero for a simple root because `mu<0` implies `mu(mu-1)>0`.

Consequently, compatibility of `q` with `A_n=R_n(O_{n-1})+S_n` is not by itself sufficient for the proposed IVT proof. It controls the sign of `A_n(mu)`, but (6) can reverse that sign.

There is a small exact counterexample satisfying the relevant degree pattern, the structural weight `C_{5,2}=4`, ordinary compatibility, and even compatibility of `q` with `A_5`; see Section 6.

---

## 4. A rigorous weighted root-separation theorem

The missing condition can be stated exactly. This theorem proves real-rootedness and interlacing with `O_{n-j}`; it does not claim consecutive interlacing.

### Notation

Assume the ordinary window

\[
 O_{n-j},O_{n-j+1},\ldots,O_{n-1} \tag{7}
\]

is commonly compatible. Set

\[
 q=O_{n-j},\qquad p=O_{n-1},\qquad e=\deg q,
 \qquad c=\binom{n-1}{j-1}.
\]

Let `mu_1<...<mu_e<0` be the roots of `q`. By Lemma 1,

\[
 \frac{p(\mu_i)}{q'(\mu_i)}<0,
 \qquad
 \frac{S_n(\mu_i)}{q'(\mu_i)}\le0. \tag{8}
\]

Define

\[
 \alpha_i=-\frac{p(\mu_i)}{q'(\mu_i)}>0,
 \qquad
 \beta_i=-\frac{S_n(\mu_i)}{q'(\mu_i)}\ge0,
 \qquad
 f_i=n+2\mu_i\frac{p'(\mu_i)}{p(\mu_i)}.
\]

Since `R_n(p)(mu_i)=p(mu_i)f_i`, equation (2) gives

\[
 \frac{O_n(\mu_i)}{q'(\mu_i)}
 =-\alpha_i f_i-\beta_i+c\mu_i(\mu_i-1). \tag{9}
\]

### Weighted root-separation inequality (WRS)

Require, for every root `mu_i` of `q`,

\[
 \boxed{
 \alpha_i f_i+\beta_i>c\mu_i(\mu_i-1).} \tag{10}
\]

Equivalently, without the auxiliary ratios,

\[
 \boxed{
 -\frac{R_n(p)(\mu_i)+S_n(\mu_i)}{q'(\mu_i)}
 >c\mu_i(\mu_i-1).} \tag{11}
\]

The right side is strictly positive.

### Theorem 2 (conditional proof of real-rootedness)

Under the common-window hypothesis (7), WRS (10) implies that `O_n` has exactly `e+1=deg O_n` distinct negative real roots. Moreover, `O_n` interlaces `q=O_{n-j}`.

**Proof.** If `e=0`, then `n=j`, WRS is vacuous, and `O_j` is linear with positive coefficients, so it has one negative root. Assume `e>=1`. By (9) and (10),

\[
 O_n(\mu_i)=-q'(\mu_i)
 \bigl(\alpha_i f_i+\beta_i-c\mu_i(\mu_i-1)\bigr),
\]

so `sgn O_n(mu_i)=-sgn q'(mu_i)`. These signs alternate at the ordered roots `mu_i`.

The degree formula is `deg O_n=floor(n/j)=e+1`, and its leading coefficient is positive. Therefore `O_n` has opposite signs at `-infinity` and `mu_1`. Consecutive values at `mu_i,mu_{i+1}` have opposite signs, and `O_n(mu_e)<0<O_n(0)=a_{n,0}`. IVT gives one root in each of

\[
 (-\infty,\mu_1),\quad
 (\mu_1,\mu_2),\ldots,\quad
 (\mu_{e-1},\mu_e),\quad
 (\mu_e,0).
\]

There are `e+1` distinct negative roots, exhausting `deg O_n`. Their positions are precisely the interlacing positions. `Q.E.D.`

Thus (10) is a genuine weighted root-preservation inequality for the exact original combination (2).

---

## 5. What the transformed-window lead does and does not prove

Suppose one proves the stronger transformed window

\[
 O_{n-j},O_{n-j+1},\ldots,O_{n-1},R_n(O_{n-1}) \tag{12}
\]

is commonly compatible. Then `A_n=R_n(O_{n-1})+S_n` is a positive sum of members other than `q`, so `q` is compatible with `A_n`. In particular,

\[
 \alpha_i f_i+\beta_i>0. \tag{13}
\]

But WRS requires the strictly larger right side `c mu_i(mu_i-1)`. Compatibility supplies only (13), not (10).

So even a complete proof of the transformed-window claim leaves a quantitative gap. A full induction must prove WRS, or a stronger invariant implying it.

For `j=2`, where `S_n=0`, WRS becomes the exact slope inequality

\[
 \boxed{
 -\frac{p(\mu)}{q'(\mu)}
 \left(n+2\mu\frac{p'(\mu)}{p(\mu)}\right)
 >(n-1)\mu(\mu-1).} \tag{14}
\]

This is the precise additional inequality that ordinary interlacing does not supply.

---

## 6. Exact obstruction: compatible inputs and compatible `A_n` are insufficient

Take `j=2`, `n=5`, so `c=4`, and choose

\[
 p(t)=(t+2)(t+4)=t^2+6t+8,
 \qquad q(t)=t+3.
\]

The roots of `p` are `-4,-2`, and the root `-3` of `q` lies strictly between them. Thus `(q,p)` is an ordinary compatible pair with exactly the relevant degrees.

Compute

\[
 A_5=R_5(p)=5p+2tp'=9t^2+42t+40.
\]

Its roots are

\[
 -\frac{10}{3},\qquad -\frac43,
\]

and `-10/3 < -3 < -4/3`. Hence `q` is also compatible with `A_5`. This example therefore satisfies both ordinary compatibility and the proposed transformed-window conclusion for `j=2`.

Nevertheless, at `mu=-3`,

\[
 p(\mu)=-1,\quad p'(\mu)=0,\quad q'(\mu)=1,
 \quad \alpha=1,\quad f=5,\quad \beta=0,
\]

while

\[
 c\mu(\mu-1)=4(-3)(-4)=48.
\]

So WRS fails: `5>48` is false. The last term is far from zero:

\[
 tQ(q)(-3)=(-3)(-4)(1)=12.
\]

The final combination is

\[
 O_5^{\rm test}=A_5+4tQ(q)
 =17t^2+50t+40,
\]

whose discriminant is

\[
 50^2-4(17)(40)=2500-2720=-220<0.
\]

Thus it has two nonreal roots.

**Conclusion.** No theorem that uses only the common compatible window, the structural coefficient `c=n-1`, and compatibility of `q` with `A_n` can prove the required root step. A quantitative invariant such as WRS is indispensable.

For comparison, the actual original-family row at `j=2,n=5` is

\[
 p=O_4=24+45t+6t^2,\qquad q=O_3=6+7t,\qquad \mu=-\frac67.
\]

Here the exact WRS margin is

\[
 -\frac{A_5(\mu)+4\mu(\mu-1)q'(\mu)}{q'(\mu)}
 =\frac{3222}{343}>0.
\]

This one exact check is illustrative only; it is not an all-index proof.

---

## 7. A refined marker recurrence, and why first moments do not close

Let `c=j-1`, `K=u\partial_u`, and mark a block of size `s<c` by `y_s`. Define

\[
 \mathcal P_n(t,u;\mathbf y)
 =\sum_{\pi}t^{b(\pi)}u^{k(\pi)}
   \prod_{s=1}^c y_s^{c_s(\pi)},
\]

where `c_s(pi)` counts blocks of size `s`. A block of size at least `j` has weight `t` and no `y` marker.

Inserting label `n` into an ordered partition of `[n-1]` gives the exact refined recurrence

\[
 \boxed{
 \begin{aligned}
 \mathcal P_n={}&K\mathcal P_{n-1}
 +u y_1(K+1)\mathcal P_{n-1}\\
 &+\sum_{s=1}^{c-1}(y_{s+1}-y_s)\partial_{y_s}\mathcal P_{n-1}
 +(t-y_c)\partial_{y_c}\mathcal P_{n-1}.
 \end{aligned}} \tag{15}
\]

**Proof.** A partition with `k` blocks has `k` choices of an existing block and `k+1` gaps for a new singleton. The existing-block choices contribute `K`, and the singleton gaps contribute `u y_1(K+1)`. Inserting into a block of size `s` changes its marker from `y_s` to `y_{s+1}`, or to `t` at the critical size. Summing that change over all blocks gives the derivative terms. `Q.E.D.`

At `u=y_1=...=y_c=1`, this recovers the original family. If

\[
 R_{n,s}=\left.y_s\partial_{y_s}\mathcal P_n\right|_{\mathbf y=1},
\]

then the specialization of (15) is

\[
 P_n(t,u)=\bigl(K+u(K+1)\bigr)P_{n-1}+(t-1)R_{n-1,c}, \tag{16}
\]

where `P_n(t,u)=sum_pi t^{b(pi)}u^{k(pi)}`. At `u=1`, the first operator becomes `1+2K`, so `K` here counts blocks, not labels.

However, the first moments `R_{n,s}` do **not** form a closed vector under (15). For example, for `j=2`,

\[
 P_1=u,\qquad R_{1,1}=u,\qquad P_2=tu+2u^2,
\]

and direct enumeration gives

\[
 R_{2,1}=4u^2.
\]

The two ordered singleton partitions each contribute two singleton blocks. The reason a closed first-moment formula fails is that inserting into a critical block changes its count from `c_s` to `c_s-1` and multiplies its weight by `t`; the correction involves factorial moments such as `sum c_s(c_s-1)t^{b(pi)}u^{k(pi)}`, or equivalently an indicator that a critical block is present. Equivalently, one must retain the full marker polynomial (15), or a larger state including presence indicators/higher marker moments.

Thus (15) is a concrete refined recurrence, but there is no positive closed first-moment production matrix: the signed threshold weight `t-1` remains, and the naive finite first-moment closure is false. No root-preserving induction for this refinement was obtained.

---

## 8. Obstruction to the proposed Simion / descent-run refinement

Simion's theorem for the total number of blocks does not automatically survive marking block sizes or marking total blocks in addition to the large-block variable.

### 8.1 Arbitrary size weights already fail

Let `P_n(x_1,...,x_n)` mark a block of size `r` by `x_r`. Its EGF is

\[
 \sum_{n\ge0}P_n\frac{z^n}{n!}
 =\frac1{1-\sum_{r\ge1}x_r z^r/r!}.
\]

For `n=4`,

\[
 P_4=x_1^4+18x_1^2x_2+6x_2^2+8x_1x_3+x_4.
\]

Specializing `x_2=x_3=x_4=1` gives

\[
 x^4+18x^2+8x+7.
\]

This polynomial has no real root. For `x>=0` it is positive. For `x=-y<0` it becomes

\[
 y^4+18y^2-8y+7>0,
\]

because `18y^2-8y+7` has discriminant `64-504=-440<0` and positive leading coefficient.

Thus even positivity of all block-size weights does not make the fully refined polynomial real-rooted. A Simion-style theorem would need a special property of the two-value specialization `(1,t)`, not merely the ordered-partition encoding.

### 8.2 Marking total blocks also fails for the two-value family

For `j=2`, mark every block additionally by `u`. The exact `n=3` polynomial is

\[
 P_3(t,u)=tu+6tu^2+6u^3.
\]

At `t=1/3`,

\[
 P_3(1/3,u)=\frac u3(18u^2+6u+1),
\]

and the quadratic factor has discriminant

\[
 6^2-4(18)(1)=-36<0.
\]

So the natural total-block refinement is not real-rooted in `u`, even though the original specialization `u=1` is the family under study. A descent-run theorem that merely marks total parts cannot be specialized naively.

No concrete Simion/multiindexed-Sturm or descent-run theorem for the required two-value size specialization was found in the supplied local sources. The exact failures above explain why such a theorem would need additional hypotheses.

---

## 9. Precise remaining target

The original recurrence (2) is valid. The corrected IVT argument reduces real-rootedness to the explicit weighted inequality WRS (10), not to unweighted compatibility.

A complete induction would need to prove, from the coefficient recurrence or from a new refined invariant, that for every `n,j` and every root `mu` of `O_{n-j}`,

\[
 -\frac{R_n(O_{n-1})(\mu)+S_n(\mu)}{O_{n-j}'(\mu)}
 >\binom{n-1}{j-1}\mu(\mu-1). \tag{19}
\]

Neither ordinary interlacing nor compatibility of `O_{n-j}` with `A_n` implies (19). Section 6 is an exact counterexample to both implications. Consecutive interlacing would require further sign inequalities at roots of the intermediate rows, beyond the sparse interlacing with `O_{n-j}` proved by Theorem 2.

## 10. Reproducible exact checks

`exact_checks.py` uses integer/rational arithmetic and exact discriminants. It checks:

- the original recurrence for small `j,n`;
- the compatible counterexample and its negative discriminant;
- the failure of WRS in that counterexample;
- the positive exact WRS margin for the actual `j=2,n=5` row;
- the refined marker recurrence for small `j,n`;
- the two refinement obstructions above.

Run from this directory:

```powershell
$py = 'C:\Users\Abdelli\AppData\Local\Python\pythoncore-3.14-64\python.exe'
& $py .\exact_checks.py
```
