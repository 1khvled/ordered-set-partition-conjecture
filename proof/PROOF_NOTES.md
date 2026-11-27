# Ordered big-block conjecture — proof-track notes

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof` was created or edited. Parent scripts were read but not changed.

## Status in one paragraph

The short recurrence guessed for `j=3` is **exactly true**, and it extends to every `j >= 2` with the same differential operator. A complete generating-function proof is given below. This is a proved identity, not a numerical guess. It does **not** prove real-rootedness or interlacing. The obvious root-compatible-window induction fails as stated: some of its component polynomials are not compatible. A second naive route — inferring real-rootedness merely from nonnegative coefficients in the positive Riccati form — is also false in general. No counterexample to the ordered big-block conjecture was found, and no complete induction or total-positivity proof was obtained.

Labels used below:

- **PROVED** — complete argument supplied here.
- **EXACT CHECK** — finite exact computation, not a proof for all `n,j`.
- **FAILED ROUTE** — a precise obstruction to the stated method.
- **OPEN** — missing argument.

---

## 1. Normalization — corrected after review

The original series already has constant term zero:
\[
 B_j(0,t)=0,\qquad F_j(0,t)=1.
\]
Indeed, the tail \(e^z-\sum_{r=0}^{j-1}z^r/r!\) vanishes at \(z=0\). Thus
\(O_n=n![z^n]F_j\) is the required integer polynomial sequence, and multiplication
of \((1-B_j)F_j=1\) by \(z^n/n!\) gives
\[
 O_0=1,\qquad
 O_n=\sum_{r=1}^{j-1}\binom nr O_{n-r}
 +t\sum_{r=j}^{n}\binom nr O_{n-r}. \tag{1}
\]

An earlier draft of this note incorrectly claimed \(B_j(0,t)=t\) and introduced a
shifted series. That claim was wrong. To avoid changing the algebra below, the
remaining sections use the harmless aliases
\[
 \widehat B_j:=B_j,\qquad G_j:=F_j.
\]
No deformation or normalization factor is being made. All recurrence and
Riccati identities below concern the original conjecture.

---

## 2. Riccati identities

### 2.1 Notation

Let

\[
 E_{j-1}(z)=\sum_{r=0}^{j-1}\frac{z^r}{r!},\qquad
 C_j(z)=e^z-E_{j-1}(z),\qquad
 w_j(z)=\frac{z^{j-1}}{(j-1)!}.
\]

Then

\[
 \widehat B_j=tC_j+(1-t)(E_{j-1}-1),
 \qquad C_j'=C_j+w_j. \tag{2}
\]

Also put

\[
 T=t\frac{\partial}{\partial t},\qquad D=\frac{\partial}{\partial t}.
\]

Since `G_j=(1-\widehat B_j)^{-1}`,

\[
 (G_j)_t=C_jG_j^2. \tag{3}
\]

A direct coefficient comparison using (2) gives

\[
 \widehat B_j'=1+\widehat B_j+(t-1)w_j. \tag{4}
\]

### 2.2 Positive form — PROVED

Differentiating `G_j=(1-\widehat B_j)^{-1}` gives

\[
 G_j'=(\widehat B_j+\widehat B_j'-1)G_j^2+G_j. \tag{5}
\]

By (4),

\[
 \widehat B_j+\widehat B_j'-1
 =2\widehat B_j+(t-1)w_j. \tag{6}
\]

The coefficients of the series in (6), in the `z/t` monomial basis, are

\[
 q_0=0,\qquad
 q_r=2\quad(1\le r\le j-2),\qquad
 q_{j-1}=1+t,\qquad
 q_r=2t\quad(r\ge j).
\]

Thus they are nonnegative for `t >= 0`. This proves the asserted positivity of the positive form, with the corrected `\widehat B_j`.

### 2.3 Short Riccati form — PROVED

Because `\widehat B_jG_j^2=G_j^2-G_j`, equations (5)–(6) imply

\[
 \boxed{G_j'=\bigl(2+(t-1)w_j\bigr)G_j^2-G_j.} \tag{7}
\]

Since `G_j=F_j` and `\widehat B_j=B_j`, this is precisely the Riccati identity for the original series `F_j=1/(1-B_j)`.

### 2.4 Coefficient form involving `H_n` — PROVED

Define

\[
 H_n=n![z^n]G_j^2
 =\sum_{s=0}^n\binom ns O_sO_{n-s}.
\]

Taking `n![z^n]` in (7) gives

\[
 \boxed{
 O_{n+1}=2H_n-O_n+(t-1)\binom n{j-1}H_{n-j+1},
 } \tag{9}
\]

where `H_k=0` for `k<0`.

---

## 3. Universal short recurrence — PROVED

### 3.1 Theorem

For every `j >= 2`, every `n >= 1`, and with negative-index terms omitted,

\[
 \boxed{
 \begin{aligned}
 O_n={}&nO_{n-1}+2T O_{n-1}
 +\sum_{r=2}^{j-1}\binom nr O_{n-r}\\
 &+\binom{n-1}{j-1}t\bigl[1+(t-1)D\bigr]O_{n-j}.
 \end{aligned}
 } \tag{10}
\]

The sum is empty for `j=2`, and the last term is absent for `n<j`.

For `j=3`, this is precisely

\[
 \boxed{
 \begin{aligned}
 O_n={}&nO_{n-1}+2tO_{n-1}'
 +\binom n2O_{n-2}\\
 &+\frac t2(n-1)(n-2)
 \bigl[1+(t-1)D\bigr]O_{n-3}.
 \end{aligned}
 } \tag{11}
\]

Thus the parent guess for `j=3` is verified exactly, including the derivative term and the binomial convention.

### 3.2 Proof of (10)

Let `G=G_j`, `C=C_j`, and `w=w_j`. Consider

\[
 \begin{aligned}
 \Phi(z,t)={}&tCG
 -2\int_0^z TG(u,t)\,du\\
 &-\frac{t}{(j-1)!}\int_0^z
 u^{j-1}\bigl[G+(t-1)G_t\bigr](u,t)\,du.
 \end{aligned} \tag{12}
\]

Clearly `\Phi(0,t)=0`. Using (2), (3), and `TG=tG_t`, differentiation gives

\[
 \partial_z\Phi
 =tC\left(G+G'-2G^2-(t-1)wG^2\right). \tag{13}
\]

The bracket vanishes by the Riccati identity (7). Hence `\Phi=0`, or

\[
 tCG
 =2\int_0^z TG(u,t)\,du
 +\frac{t}{(j-1)!}\int_0^z
 u^{j-1}\bigl[G+(t-1)G_t\bigr](u,t)\,du. \tag{14}
\]

Now extract `n![z^n]`. The left side is

\[
 t\sum_{r=j}^n\binom nr O_{n-r}.
\]

The first integral contributes `2T O_{n-1}`. For the second integral, put `m=n-j`; then

\[
 n![z^n]\int_0^z u^{j-1}P(D)G(u,t)\,du
 =\frac{(n-1)!}{(n-j)!}P(D)O_{n-j},
\]

and therefore

\[
 t\sum_{r=j}^n\binom nr O_{n-r}
 =2T O_{n-1}
 +\binom{n-1}{j-1}t\bigl[1+(t-1)D\bigr]O_{n-j}. \tag{15}
\]

Finally substitute (15) into (1) and separate the `r=1` term. This is exactly (10). `Q.E.D.`

### 3.3 Coefficient form useful for PF/TP work — PROVED

Write `a_{n,m}=[t^m]O_n`, with out-of-range entries zero. Since

\[
 [t^m]\,t\bigl[1+(t-1)D\bigr]O
 =m(a_{m-1}-a_m),
\]

equation (10) is equivalent to

\[
 \boxed{
 \begin{aligned}
 a_{n,m}={}&(n+2m)a_{n-1,m}
 +\sum_{r=2}^{j-1}\binom nr a_{n-r,m}\\
 &+\binom{n-1}{j-1}m
 \bigl(a_{n-j,m-1}-a_{n-j,m}\bigr).
 \end{aligned}
 } \tag{16}
\]

The negative term in (16) is one reason a routine positive-production-matrix proof is not immediate.

### 3.4 Degrees and coefficient signs — PROVED

For `n >= 1`,

\[
 \deg_t O_n=\left\lfloor\frac nj\right\rfloor,
 \qquad a_{n,m}>0\quad
 \left(0\le m\le\left\lfloor\frac nj\right\rfloor\right). \tag{17}
\]

Indeed, `m` big blocks require at least `mj` labels, giving the degree bound. Conversely, for any `m <= floor(n/j)`, choose `m` disjoint blocks of size `j` and make all remaining labels singleton blocks; ordering the blocks gives a positive contribution. Since all coefficients are positive, every real root of `O_n` is negative.

---

## 4. Exact computational verification

All arithmetic in the verification scripts is integer or rational arithmetic. No floating-point root decision is used for the identities.

`verify_identities.py` checks, for `j=2,...,12` and `n<=60`:

- coefficient generation from (1);
- the short recurrence (10);
- the `H_n` Riccati coefficient identity (9);
- the degree formula (17).

Every check passes. Results are in `exact_results.json`.

These finite checks are consistency tests for the proved identities, not the proof of the root conjecture.

---

## 5. Root-preserving pieces of a possible induction

### 5.1 A useful operator — PROVED

Let `p` be real-rooted with positive leading coefficient and only negative roots, of degree `d`. For `n>0`, define

\[
 R_n(p)=np+2t p'.
\]

Then `R_n(p)` has `d` negative real roots. More precisely, if

\[
 r_1<\cdots<r_d<0
\]

are the roots of `p`, the roots `s_i` of `R_n(p)` satisfy

\[
 r_i<s_i<r_{i+1}\quad(i<d),
 \qquad r_d<s_d<0. \tag{18}
\]

**Proof.** Write `p(t)=c\prod_{k=1}^d(t+\alpha_k)`, `\alpha_k>0`. Away from the roots of `p`, a root of `R_n(p)` is a zero of

\[
 f(t)=n+2\frac{t p'(t)}{p(t)}
 =n+2\sum_{k=1}^d\frac{t}{t+\alpha_k}.
\]

On every interval between consecutive roots, each summand has derivative

\[
 \frac{\alpha_k}{(t+\alpha_k)^2}>0.
\]

At the left endpoint of such an interval, the term belonging to that endpoint tends to `-\infty`; at the right endpoint it tends to `+\infty`. Thus there is exactly one zero in each gap. On `(r_d,0)`, `f` goes from `-\infty` to `n>0`, giving one final zero. These are `d` distinct roots, which exhaust the degree. `Q.E.D.`

Thus the first two terms of (10) individually preserve negative real-rootedness and move the roots to the right.

### 5.2 The last differential operator — PROVED, but with a caveat

Define

\[
 Q(p)=p+(t-1)p'=\bigl((t-1)p\bigr)'.
\]

If `p` has `d` negative real roots, then `(t-1)p` has those roots together with the root `1`. Rolle's theorem therefore shows that `Q(p)` has `d` real roots, one in each gap.

However, the largest root lies in `(r_d,1)` and **need not be negative**. So `Q` preserves real-rootedness but not negative-rootedness.

### 5.3 Positive combinations of compatible polynomials — standard proved closure

If two real-rooted polynomials have adjacent degrees and their roots interlace in the appropriate orientation, then every positive linear combination is real-rooted. The same holds for equal-degree polynomials whose roots strictly alternate with the correct endpoint order. This is the classical Obreschkoff/compatible-polynomial closure; it follows by applying Rolle's theorem to the Wronskian, or equivalently by observing that `-p/q` is strictly monotone on the relevant root intervals.

This gives a valid conditional induction **if** all terms in (10) can be placed in one compatible window with positive coefficients.

---

## 6. FAILED ROUTE: the literal compatible window is not compatible

Take `j=2`, `n=5`. The needed polynomials are

\[
 O_3=6+7t,\qquad O_4=24+45t+6t^2.
\]

The two differential terms in (10) are

\[
 R_5(O_4)=5O_4+2TO_4
 =120+315t+54t^2,
\]

and

\[
 Q(O_3)=O_3+(t-1)O_3'
 =14t-1.
\]

The polynomial `O_4` has two negative roots because its discriminant is

\[
 45^2-4(6)(24)=1449>0
\]

and both the sum and product of its roots are positive. But `Q(O_3)` has the positive root `1/14`.

Since `deg O_4=2` and `deg Q(O_3)=1`, compatibility would require the root of `Q(O_3)` to lie between the two negative roots of `O_4`. It does not. The pair `(R_5(O_4),Q(O_3))` is likewise incompatible: `R_5(O_4)` also has two negative roots.

`window_probe.py` confirms this with rigorous Arb enclosures and finds the same obstruction repeatedly for small `j,n`; the first exact witness is the one above.

This is a genuine obstruction to the **unweighted pairwise-compatible-window method**, not to the conjecture. Indeed the actual weighted sum is

\[
 O_5=R_5(O_4)+4tQ(O_3)
 =120+311t+110t^2,
\]

whose discriminant is `43921>0` and whose roots are negative.

A successful induction would need a weighted compatibility theorem that uses the actual coefficients in (10), rather than mere pairwise compatibility of all terms.

---

## 7. FAILED ROUTE: nonnegative Riccati coefficients alone are insufficient

The positive form (5) has nonnegative `z/t` coefficients for the actual family. It is tempting to infer real-rootedness from that positivity alone. The implication is false.

Consider instead

\[
 G'=t^2zG^2+G,\qquad G(0,t)=1.
\]

The potential `t^2z` has nonnegative coefficients. If

\[
 O_n=n![z^n]G,
\]

then

\[
 O_0=1,\qquad O_1=1,\qquad O_2=1+t^2.
\]

The last polynomial has discriminant `-4` and no real root. This exact counterexample is verified by `riccati_counterexample.py`.

Therefore any proof must use additional structure of the actual potential `2\widehat B_j+(t-1)w_j`, not positivity by itself.

---

## 8. PF/TP route: evidence and missing structure

Let

\[
 A^{(j)}_{n,m}=a_{n,m}.
\]

A standard sufficient strategy is to prove that this graded coefficient matrix is totally positive, then apply the interlacing-polynomial/total-positivity criterion. The degree staircase (17) must be respected when forming minors.

`minor_probe.py` performed a bounded exact check of all minors of sizes `2` through `5` for `j=2,...,6` and rows `0,...,18`. No negative minor was found. This is evidence only.

What is still missing for a proof is one of the following:

1. an explicit oscillatory or sign-regular production matrix for `A^{(j)}`;
2. a planar-network/Lindström–Gessel–Viennot interpretation of `a_{n,m}` giving nonintersecting-path determinants;
3. a direct proof that the relevant minors satisfy a positive recurrence;
4. a replacement for total positivity tailored to the weighted window in (10).

Equation (16) does not by itself supply a positive banded production matrix: it contains the difference `a_{n-j,m-1}-a_{n-j,m}`. No valid infinite-minor argument was obtained in this track.

---

## 9. Audit warning about the parent interlacing checker

This is a read-only observation about `ordered-big-blocks/interlace.py`; it was not edited.

In `check_pair`, `a` and `b` are lists of roots, so their lengths are the polynomial degrees. The code then sets

```python
d = len(a) - 1
e = len(b) - 1
```

For the degree-jump branch this off-by-one happens to cancel with the loop range. The equal-degree branch, however, checks only

```python
all(a[i] <= b[i] and b[i] <= a[i + 1] for i in range(d))
```

and omits the final root of `b` and the final endpoint inequality. For linear polynomials `d=0`, so the loop performs no check at all. Thus equal-degree interlacing can be falsely certified.

A correct equal-degree check must compare all `len(a)` roots and include the final inequality in the chosen orientation. This does not indicate a counterexample to the conjecture, but it means equal-degree entries in the parent audit should not be treated as certificates without fixing or independently checking that branch.

---

## 10. Precisely what remains open

The universal recurrence (10) and all Riccati identities in this note are proved. The following are not:

- a base-to-infinity proof that every `O_{n,j}` has only negative real roots;
- a proof of consecutive interlacing, especially across degree jumps at multiples of `j`;
- a weighted replacement for the failed compatible-window argument;
- a total-positivity or PF production-matrix proof;
- a counterexample to the actual ordered big-block conjecture.

The most concrete next targets are:

1. use (16) to derive a signed recurrence for the `2x2` PF minors and prove their positivity;
2. identify a weighted compatibility inequality for the exact combination in (10);
3. construct an explicit planar path model for `a_{n,m}`;
4. repair the equal-degree branch of the parent interlacing audit before relying on its output.

Until one of these is completed, the root conjecture remains open in this proof track.

---

## 11. Files in this proof directory

- `PROOF_NOTES.md` — this note.
- `exact_algebra.py` — integer-polynomial arithmetic, coefficient generation, and recurrence residual.
- `verify_identities.py` — bounded exact checks; writes `exact_results.json`.
- `window_probe.py` — rigorous root-enclosure compatibility probe; writes `window_results.json`.
- `minor_probe.py` — bounded exact coefficient-minor probe; writes `minor_results.json`.
- `riccati_counterexample.py` — exact obstruction to the positivity-only Riccati route.
- `exact_results.json`, `window_results.json`, `minor_results.json` — reproducible finite outputs.

Reproduction commands (run from the proof directory):

```powershell
$py = 'C:\Users\Abdelli\AppData\Local\Python\pythoncore-3.14-64\python.exe'
& $py .\verify_identities.py
& $py .\window_probe.py
& $py .\minor_probe.py
& $py .\riccati_counterexample.py
```
