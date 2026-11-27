# Literature and theory memo: ordered large-block polynomials

**Checked:** 27 November 2026  
**Family:** for fixed $j\ge2$,
\[
O_{n,j}(t)=\sum_{\pi\text{ ordered partition of }[n]}
 t^{\#\{B\in\pi:|B|\ge j\}},
\qquad
F_j(z,t)=\sum_{n\ge0}O_{n,j}(t)\frac{z^n}{n!}
 =\frac1{1-B_j(z,t)},
\]
where
\[
B_j(z,t)=A_j(z)+tC_j(z),\qquad
A_j(z)=\sum_{r=1}^{j-1}\frac{z^r}{r!},\qquad
C_j(z)=e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}.
\]

## Executive conclusion

1. **The exact enumerative family is already in the literature.** Chen and Zhu's 2026 paper *A Snail Race Problem* defines $r_m(n,k)$ exactly as these coefficients. Its Theorem 3 is the correct first-block recurrence with $\binom ni$, and Theorem 7 is the same EGF. That paper does not prove real-rootedness, interlacing, or total positivity.
2. **The all-$j$ real-rootedness/consecutive-interlacing conjecture is still open.** It is stated by Alexandersson on the current SymCat open-problem page and tracked as [solved-conjectures#108](https://github.com/1337isnot1337/solved-conjectures/issues/108). The repository now contains partial edge-zero theorems, not a complete proof.
3. **Alexandersson's proved “ordered lift” theorem is a different family.** It orders only the large blocks of partitions whose other blocks are singletons, with EGF
   \[
   \frac{e^z}{1-tC_k(z)}.
   \]
   Ordering **all** blocks gives
   \[
   \frac1{1-A_j(z)-tC_j(z)}.
   \]
   The two EGFs, and the two combinatorial families, are not equal.
4. **The coefficient matrix is not literally strictly totally positive**, because it has structural zeros. The natural matrix statement is total nonnegativity ($\mathrm{TN}_\infty$: every minor $\ge0$), or the stronger support-strict version in which every minor not forced to vanish by the support is $>0$.
5. **No published theorem proving either matrix statement for the target family was located.** Standard proper-Riordan and Pólya-frequency criteria do not apply directly. The most relevant general left-production theorem has an explicit negative-minor obstruction in this case, so that proof route cannot simply be invoked.
6. **Exact finite evidence is strong but finite.** The accompanying checker verifies every minor of sizes $2$ through $10$ in the leading $10\times10$ coefficient matrix for $j=2,3,4$: $184655$ minors per $j$, all nonnegative. Every support-admissible minor in those checks is strictly positive. This is not an all-index proof.

---

## 1. The correct recurrence and coefficient formula

Set $O_{0,j}=1$. In an ordered partition, choose its **first block**. If that block has size $r$, there are $\binom nr$ ways to choose its labels. It contributes weight $1$ for $r<j$ and weight $t$ for $r\ge j$. Hence
\[
\boxed{
O_{n,j}(t)=
\sum_{r=1}^{j-1}\binom nr O_{n-r,j}(t)
+t\sum_{r=j}^{n}\binom nr O_{n-r,j}(t)
}
\tag{1}
\]
with sums restricted to $r\le n$ and negative-index rows zero.

The factor is $\binom nr=\binom n{n-r}$, **not** $\binom{n-1}{r-1}$. The latter is the insertion recurrence for an *unordered* set partition with a distinguished new label; it does not choose the first block of an ordered partition from all $n$ labels. Already for $j=2,n=2$, the correct recurrence gives
\[
O_{2,2}=2+t,
\]
whereas the wrong binomial choice gives $1+t$ and misses the second ordering of the two singleton blocks.

Writing $a_{n,m}=[t^m]O_{n,j}(t)$, geometric expansion in $tC_j/(1-A_j)$ gives
\[
\boxed{
a_{n,m}=n![z^n]\frac{C_j(z)^m}{(1-A_j(z))^{m+1}}.
}
\tag{2}
\]
In particular,
\[
a_{n,m}=0\quad(m\ge1,\ n<jm),\qquad
\deg O_{n,j}=\left\lfloor\frac nj\right\rfloor,
\]
and $a_{n,0}>0$ for every $n\ge0$.

The issue #108 research notes also derive
\[
F_j'=\left(2+(t-1)\frac{z^{j-1}}{(j-1)!}\right)F_j^2-F_j.
\tag{3}
\]
If
\[
H_n=n![z^n]F_j^2
 =\sum_{r=0}^n\binom nr O_{r,j}O_{n-r,j},
\]
then (3) is equivalent to
\[
O_{n+1,j}=2H_n-O_{n,j}
 +(t-1)\binom n{j-1}H_{n-j+1,j}.
\tag{4}
\]
The signed $(t-1)$ term is the obstruction to a direct positive-convolution/root-preservation induction.

---

## 2. Direct prior literature for the exact coefficients

### Chen–Zhu, *A Snail Race Problem* (2026)

**Source:** Yu Chen and Wenxiang Zhu, [arXiv:2607.02085](https://arxiv.org/html/2607.02085).

The paper defines $r_m(n,k)$ as the number of snail-race outcomes with exactly $k$ ties of size at least $m$. An outcome is an ordered partition, so
\[
r_m(n,k)=a_{n,k}\quad\text{for }m=j.
\]

**Theorem 3 (exact statement).** For $m,n\in\mathbb N$ and $k\ge0$,
\[
r_m(n,k)=
\sum_{i=1}^{m-1}\binom ni r_m(n-i,k)
+\sum_{i=m}^{n}\binom ni r_m(n-i,k-1),
\]
with $r_m(0,0)=1$ and the stated zero boundary conventions. This is exactly (1).

**Theorem 7 (exact statement).** Let
\[
T_{m-1}(x)=\sum_{q=0}^{m-1}\frac{x^q}{q!}.
\]
Then the EGF in $n$ of $r_m(n,k)$ is
\[
\sum_{n\ge0}r_m(n,k)\frac{x^n}{n!}
=
\frac{[e^x-T_{m-1}(x)]^k}
     {[2-T_{m-1}(x)]^{k+1}}.
\]
Since $2-T_{m-1}=1-A_j$ and $e^x-T_{m-1}=C_j$, this is exactly (2).

The paper also proves closed formulas and connects the numbers to ordered Bell, Stirling, and partial Bell polynomials. **It contains no zero theorem and no total-positivity theorem.**

For $j=2,k=1$, the resulting column is
\[
1,7,45,311,2383,20301,\ldots
\]
(from $n=2$ onward), recorded as [OEIS A397265](https://oeis.org/A397265). The OEIS entry gives the same correct recurrence and cites Chen–Zhu; it does not state a root or matrix-positivity result.

---

## 3. Ordinary ordered Bell/Fubini polynomials: the proved baseline

If every block is counted, i.e. the threshold is $j=1$, then
\[
F(z,t)=\frac1{1-t(e^z-1)},\qquad
O_{n,1}(t)=\sum_{k=0}^n k!\,S(n,k)t^k.
\tag{5}
\]
These are the ordinary ordered Bell/Fubini polynomials.

**Known facts.**

- Comtet's recurrence is
  \[
  O_{n+1,1}=t(1+D)O_{n,1},\qquad O_{1,1}=t.
  \]
  It gives only real zeros in $[-1,0]$. Consecutive rows interlace in the usual weak sense; the common zero at $0$ disappears after dividing the positive-degree rows by $t$.
- The coefficient triangle $k!S(n,k)$ is a positive column scaling of the Stirling triangle $S(n,k)$. The Stirling triangle is totally nonnegative by Chen–Fu–Ruan, Corollary 4.3, so the Fubini triangle is $\mathrm{TN}_\infty$.
- Jalowy, Kabluchko, and Marynych, *Zeros and Exponential Profiles of Polynomials II: Examples*, [arXiv:2509.11248](https://arxiv.org/html/2509.11248), published as [Constructive Approximation (2026), DOI 10.1007/s00365-026-09775-2](https://link.springer.com/article/10.1007/s00365-026-09775-2), record the classical real-rootedness in §4.3 and prove in Theorem 4.4 the limiting empirical zero distribution
  \[
  \mu_F\quad\text{on }(-1,0),\qquad
  G_F(t)=\frac1{t(t+1)\log(1+1/t)}.
  \]
  This is current asymptotic literature for the **all-blocks-counted** family, not for $O_{n,j}$ with $j\ge2$.

Thus the ordinary Fubini case is fully understood at the coefficient-matrix level. The difficulty begins when small blocks are left uncounted.

---

## 4. The ordered-lift theorem that must not be conflated with the target

Fix $k\ge2$ and let $S_{n,k}(t)$ count **unordered** set partitions of $[n]$ whose blocks have size $1$ or at least $k$, with $t$ counting the large blocks. Then
\[
\sum_{n\ge0}S_{n,k}(t)\frac{z^n}{n!}
=\exp\!\left(z+tC_k(z)\right).
\tag{6}
\]
Write $S_{n,k}(t)=\sum_r s_{n,k,r}t^r$. Alexandersson's ordered lift is
\[
\widehat S_{n,k}(t)=\sum_r r!\,s_{n,k,r}t^r.
\tag{7}
\]
It orders only the $r$ large blocks; the singleton blocks remain an unordered set. Therefore
\[
\boxed{
\sum_{n\ge0}\widehat S_{n,k}(t)\frac{z^n}{n!}
=\frac{e^z}{1-tC_k(z)}.
}
\tag{8}
\]

**Theorem (Alexandersson, 2026), as stated on SymCat.** For every $k\ge2$ and $n\ge0$, both $\widehat S_{n,k}$ and $S_{n,k}$ have only real nonpositive zeros, and
\[
\widehat S_{n,k}\interl\widehat S_{n+1,k},
\qquad
S_{n,k}\interl S_{n+1,k}.
\]
The stated proof uses a sliding window consisting of a row and its Euler transform, preserved by Wagner's interlacing cones; the multiplier sequence $(1/r!)$ then removes the ordering.

**Source:** Per Alexandersson, [Real-rooted polynomials from set partitions, “Big blocks”](https://www.symmetricfunctions.com/realRootedSetPartitions.htm#bigBlockSetPartitions), page dated 18 November 2026. This is an author-maintained theorem statement; the page does not give a separate peer-reviewed preprint for this specific theorem.

### Why (8) is not the target

The target orders **every** block, including the small ones, and allows all small sizes $1,\ldots,j-1$. Its EGF is
\[
\boxed{
F_j(z,t)=\frac1{1-A_j(z)-tC_j(z)}
=\frac1{1-A_j(z)}
 \frac1{1-tC_j(z)/(1-A_j(z))}.
}
\tag{9}
\]
The factor $1/(1-A_j)$ in (9) is an **ordered sequence** of small blocks. It accounts for arbitrary strings of small blocks before, between, and after large blocks. In (8), $e^z$ is an **unordered set** of singleton blocks, and only the large blocks are ordered.

Even for $k=j=2$,
\[
\widehat F_{2}(z,t)=\frac{e^z}{1-t(e^z-1-z)},
\qquad
F_{2}(z,t)=\frac1{1-z-t(e^z-1-z)},
\]
and already $O_{2,2}=2+t$ while the ordered lift has constant term $1$ at $n=2$. For $j\ge3$, the target additionally orders blocks of sizes $2,\ldots,j-1$. Consequently, Alexandersson's theorem proves a neighboring family, not $O_{n,j}$.

---

## 5. Current status of the target conjecture and issue #108

**Primary statement:** [Alexandersson's open-problem page, “Set partitions and trees”](https://www.symmetricfunctions.com/realRootedOpenProblems.htm#set-partitions-and-trees), checked 27 November 2026. It states the all-$j$ conjecture for $O_{n,j}$ and records finite evidence for $j\le5$, $n\le30$.

**Issue:** [1337isnot1337/solved-conjectures#108](https://github.com/1337isnot1337/solved-conjectures/issues/108), still open. The cited commits were inspected read-only:

- [`04de706`](https://github.com/1337isnot1337/solved-conjectures/commit/04de7069df504ff494d81085f3bbb6b0d62bb3ed): the correct first-block recurrence, coefficient formula (2), Riccati identity (3), and convolution recurrence (4).
- [`e76826a`](https://github.com/1337isnot1337/solved-conjectures/commit/e76826ae91fb6415083377ee25c5505e6462b272): correction of the Newton-inequality normalization in the finite diagnostic code.
- [`f1a5ac2`](https://github.com/1337isnot1337/solved-conjectures/blob/f1a5ac298febcfd20b19e1ef3465325eb62e18d6/research/alexandersson-ordered-set-partitions/edge.md): a dominant-pole argument claiming simple negative zeros and strict interlacing for a positive fraction of the rightmost real zeros, eventually for each fixed $j$.
- [`e39bd01`](https://github.com/1337isnot1337/solved-conjectures/blob/e39bd018c0f980ff394daa7764272226890c6d14/research/alexandersson-ordered-set-partitions/left-edge.md): a reciprocal-polynomial limit claiming eventual simple negative zeros and cyclic interlacing for each fixed number of zeros of largest absolute value.

These are repository research proofs, not a published solution of the full conjecture. They leave the middle zeros and the all-row statement open.

A second current tracking source is [PerAlexandersson/RealRooted#2590](https://github.com/PerAlexandersson/RealRooted/issues/2590), which explicitly lists “the ordered big-block family $1/(1-B_j)$” as open. Its newer exact two-block-size classification ([#2665](https://github.com/PerAlexandersson/RealRooted/issues/2665)) does not apply: the target has unmarked sizes $1,\ldots,j-1$ and an infinite marked tail. Even $j=2$ has one unmarked size but infinitely many marked sizes, not two block sizes.

---

## 6. What total positivity would mean here

There are two conventions in the literature:

- Chen, Zhu, and Chen–Fu–Ruan use **TP** to mean all minors are nonnegative.
- Brändén, SymCat, and much of matrix analysis use **TN** for nonnegative minors and reserve **TP** for strictly positive minors.

This memo uses **TN** for $\ge0$ and **support-strict TN** for “$>0$ whenever the minor is not forced to be zero by the support.”

The support of the target matrix is
\[
a_{n,m}>0
\iff
m=0\ \text{or}\ n\ge jm.
\tag{10}
\]
Thus literal strict TP is impossible. The meaningful conjectural matrix statement is
\[
\det a_{I,J}\ge0
\quad\text{for all finite increasing }I,J,
\tag{11}
\]
with the stronger support-strict version as a possible proof target.

A caution: **TN of the coefficient array alone does not imply real-rooted rows.** For example, the lower-triangular all-ones staircase matrix is TN, but its row $1+x+x^2$ is not real-rooted. Row Pólya frequency or a suitable production/network theorem must also be proved.

---

## 7. Riordan and production-matrix theorems: exact statements and applicability

### 7.1 Xi Chen–Liang–Wang: proper Riordan arrays

**Source:** [*Total positivity of Riordan arrays*, arXiv:1601.05637](https://arxiv.org/html/1601.05637), European Journal of Combinatorics 46 (2015), 68–74.

**Theorem 2.1.** Let $R=(r_{n,k})$ be the proper Riordan array $(g,f)$ with $g(0)=1$, $f(0)=0$, $f'(0)\ne0$, and let
\[
J(R)=
\begin{pmatrix}
z_0&a_0&&\\
z_1&a_1&a_0&\\
z_2&a_2&a_1&a_0\\
\vdots&\vdots&\vdots&\ddots
\end{pmatrix}
\]
be its coefficient/production matrix. If $J(R)$ is $\mathrm{TP}_r$ (TP), then $R$ is $\mathrm{TP}_r$ (TP).

**Applicability gap.** Remove the harmless row scaling $a_{n,m}=n!b_{n,m}$ and put
\[
d(z)=\frac1{1-A_j(z)},\qquad f(z)=C_j(z)d(z).
\]
Then
\[
b_{n,m}=[z^n]d(z)f(z)^m.
\]
But $f$ has order $j$, so $f'(0)=0$ for every $j\ge2$. This is a **degenerate**, not proper, Riordan array. Theorem 2.1 cannot be invoked without a new extension to degenerate arrays.

### 7.2 Bao-Xuan Zhu: exponential Riordan arrays

**Source:** [*Total positivity from the exponential Riordan arrays*, arXiv:2006.14485v3](https://arxiv.org/html/2006.14485v3), SIAM Journal on Discrete Mathematics 35 (2021), 2971–3003.

**Theorem 3.6.** Suppose $g$ is a $\mathrm{PF}_{r_1}$ function. If either $f$ or $1/\bar f\,'$ is $\mathrm{PF}_{r_2}$, then
\[
(g e^{\lambda f},f),\qquad (g(f)e^{\lambda f},f)
\]
are $\lambda$-$\mathrm{TP}_r$, where $r=\min(r_1,r_2)$ and $\bar f$ is the compositional inverse of $f$.

**Applicability gaps.**

1. The definition requires $f_1=f'(0)\ne0$; the target has $f'(0)=0$.
2. The needed PF hypotheses are not available for the natural factorization. For $j\ge3$, $d=1/(1-A_j)$ is not even $\mathrm{PF}_2$: its first coefficients are $1,1,3/2$, and the Toeplitz $2\times2$ minor on rows $\{0,1\}$ and columns $\{1,2\}$ is
   \[
   1^2-1\cdot\frac32=-\frac12.
   \]
3. The tail kernel also fails PF; see §8.2.

### 7.3 Chen–Fu–Ruan: left production matrices, including noninvertible arrays

**Source:** Xi Chen, Lang Fu, and Jiajie Ruan, [*A unified planar network approach to total positivity of combinatorial matrices and real-rootedness of polynomials*, arXiv:2512.08369](https://arxiv.org/html/2512.08369), December 2025.

**Theorem 1.3.** Let $A=(a_{n,k})$ be lower triangular, possibly noninvertible, and let
\[
A=Q(A)\operatorname{diag}(I_1,Q(A))
       \operatorname{diag}(I_2,Q(A))\cdots
\]
be its left-production factorization. If $Q(A)$ is TN, then:

1. $A$ is TN;
2. the reversal $A^{\leftarrow}_{n,k}=a_{n,n-k}$ is TN;
3. every row-generating polynomial $\sum_k a_{n,k}x^k$ is real-rooted.

**Theorem 1.4.** The leading Toeplitz matrix of row $n$ is the transpose of an explicitly described submatrix of a product of block copies of the leading truncation of $Q(A)$. This is the mechanism proving row Pólya frequency in Theorem 1.3.

**Corollary used in the RealRooted formalization.** If
\[
F(z,x)=\frac{d(z)}{1-xz\,d(z)K(z)},
\]
and the coefficient sequences of both $d$ and $K$ are PF, then the coefficient array is TN and its row polynomials have only real nonpositive roots. The left production matrix is
\[
T(d)\operatorname{diag}(1,T(K)).
\]

This is the closest general theorem to the target because it allows noninvertible arrays. It still does **not** apply through its standard PF corollary, and its left-production hypothesis is actually false for the target's natural factorization; see §8.

Also note: Theorem 1.3 proves **real-rootedness**, not consecutive interlacing. Interlacing requires an additional strict-TP, Wronskian, or root-sign argument.

---

## 8. Explicit obstruction to the standard PF/left-production route

### 8.1 The natural degenerate Riordan factorization

Define
\[
K_j(z)=\frac{C_j(z)}z
=\sum_{q\ge j-1}\frac{z^q}{(q+1)!}.
\]
Then
\[
b_{n,m}=[z^n]d(z)^{m+1}(zK_j(z))^m,
\qquad d(z)=\frac1{1-A_j(z)}.
\tag{12}
\]
The Chen–Fu–Ruan left production matrix is therefore
\[
Q_{n,0}=d_n,\qquad
Q_{n,m}=[z^{n-m}]d(z)K_j(z)\quad(m\ge1).
\tag{13}
\]

### 8.2 Neither natural PF hypothesis holds

For every $j\ge2$, $K_j$ is not PF. A exact negative Toeplitz minor is
\[
\boxed{
\det T(K_j)[\{j,j+1,j+2\},\{0,1,2\}]
=
\frac{1-j}{(j!)^3(j+1)^3(j+2)(j+3)}<0.
}
\tag{14}
\]
For $j=2$ this is $-1/4320$.

For every $j\ge3$, $d_j=1/(1-A_j)$ is not even $\mathrm{PF}_2$, as shown by the minor $-1/2$ above. For $j=2$, $d_2=1/(1-z)$ is PF, but $K_2$ still fails (14).

Therefore the PF-Riordan corollary cannot be used for any $j\ge2$.

### 8.3 The left production matrix itself is not TN

The failure is stronger than a missing PF certificate. For $j=2$, rows $\{2,3,4,5\}$ and columns $\{0,1,2,3\}$ of $Q$ give
\[
\begin{pmatrix}
1&1/2&0&0\\
1&2/3&1/2&0\\
1&17/24&2/3&1/2\\
1&43/60&17/24&2/3
\end{pmatrix},
\qquad
\det=-\frac1{4320}.
\tag{15}
\]
For $j=3$, rows $\{3,4\}$ and columns $\{0,1\}$ give
\[
\begin{pmatrix}
2&1/6\\
11/4&5/24
\end{pmatrix},
\qquad
\det=-\frac1{24}.
\tag{16}
\]

Thus the hypothesis “$Q(A)$ is TN” in Chen–Fu–Ruan Theorem 1.3 fails for the natural factorization. **This is not a counterexample to (11)**: it only blocks that particular sufficient theorem. A different network, factorization, or direct minor proof could still establish TN of the target matrix.

---

## 9. Reciprocal bivariate series and Pólya-frequency kernels

### 9.1 Brändén–Saud Maia Leite: the nearest reciprocal-series theorems

**Source:** [*Totally nonnegative matrices, chain enumeration and zeros of polynomials*, arXiv:2412.06595v3](https://arxiv.org/html/2412.06595v3), December 2025.

**Theorem 4.4.** Let $f$ be a Pólya-frequency series in the Aissen–Edrei–Schoenberg–Whitney form, and define
\[
\frac1{1-t(f(x)-f(0))}=\sum_{n\ge0}r_n(t)x^n.
\]
Then every $r_n$ is real-rooted and $r_n\interl r_{n+1}$. If $f(0)\ne0$, all zeros lie in $[-1/f(0),0]$.

**Theorem 4.5 (Forgács–Tran problem).** Let $Q$ be a polynomial with $Q(0)>0$ whose zeros are all real and positive, and let $r\ge1$. Define
\[
\frac1{Q(x)-t x^r}=\sum_{n\ge0}q_n(t)x^n.
\]
Then every $q_n$ has only real negative zeros and $q_n\interl q_{n+1}$.

**Why these do not solve the target.**

- Theorem 4.4 requires the whole nonconstant part of the denominator to be multiplied by the same $t$. The target has an uncounted term $A_j(z)$ independent of $t$.
- Theorem 4.5 requires the marked term to be a monomial $t x^r$. Here it is $tC_j(z)$, an infinite tail.
- For $j\ge3$, $1-A_j(z)$ does not have only positive real zeros. For example, $1-z-z^2/2$ has a negative zero.
- The natural kernel $C_j(z)/z=K_j(z)$ is not PF by (14).

The Ahlberg–Alexandersson weighted-renewal theorem on [SymCat's interlacing page](https://www.symmetricfunctions.com/realRootedInterlacing.htm#interlacingRenewalRecurrences) is also different: it treats
\[
\frac1{1-tR(xt^\ell)}
\]
with a real-rooted polynomial seed $R$. The target has $C_j(z)/(1-A_j(z))$ in place of $R(xt^\ell)$, and that kernel is neither a real-rooted polynomial seed nor PF.

### 9.2 Pólya-frequency sequences versus continuous Pólya-frequency kernels

**Finite/discrete PF criterion.** A polynomial $\sum_{k=0}^d a_k t^k$ with nonnegative coefficients has only real nonpositive zeros iff its coefficient sequence, extended by zeros, is a PF sequence, i.e. its lower Toeplitz matrix is TN.

**Source:** Aissen, Edrei, Schoenberg, and Whitney, [*On the Generating Functions of Totally Positive Sequences*, PNAS 37 (1951), 303–307](https://doi.org/10.1073/pnas.37.5.303), especially Theorems 4 and 6.

**Continuous PF-kernel theory.** A Pólya-frequency function $\Lambda$ is an integrable function, nonzero at two or more points, such that the Toeplitz kernel
\[
T_\Lambda(x,y)=\Lambda(x-y)
\]
is TN. Schoenberg's theorem identifies its bilateral Laplace transform with the reciprocal of a function in the Laguerre–Pólya class.

**Source:** Belton, Guillot, Khare, and Putinar, [*Totally positive kernels, Polya frequency functions, and their transforms*, arXiv:2006.16213](https://arxiv.org/abs/2006.16213), Journal d'Analyse Mathématique 150 (2023), 83–158; see Definition 8.1 and the discussion of Schoenberg's transform theorem.

This continuous theory is useful background for variation diminution and reciprocal Laplace transforms. It is **not** a theorem that an arbitrary formal reciprocal bivariate series has a TN coefficient array. No analytic continuation/separability hypothesis needed to apply it has been established for $F_j(z,t)$.

---

## 10. Is the target coefficient matrix known or provably TN?

### 10.1 Literature status

No source located in this search states or proves
\[
\det(a_{n_i,m_j})_{i,j=1}^r\ge0
\]
for all finite increasing row and column sets, for any fixed $j\ge2$. In particular:

- Chen–Zhu prove the enumeration, not positivity.
- Alexandersson's ordered-lift theorem concerns a different EGF.
- The proper Riordan criteria exclude the target because $f'(0)=0$.
- The general Chen–Fu–Ruan criterion is blocked by the explicit negative minors (15)–(16) in its left production matrix.
- The Brändén–Leite reciprocal theorems require PF or monomial/positive-root hypotheses that fail here.

### 10.2 Exact finite evidence from this memo

The accompanying script is
[`literature/check_coefficient_minors.js`](literature/check_coefficient_minors.js).
Run from `ordered-big-blocks` with:

```text
node literature/check_coefficient_minors.js
```

It generates rows by (1) using exact `BigInt` binomial coefficients and evaluates minors by fraction-free Bareiss elimination. Output:

```json
{"j":2,"side":10,"tested":184655,"admissible":1397,"negative":null,"admissibleZero":null}
{"j":3,"side":10,"tested":184655,"admissible":317,"negative":null,"admissibleZero":null}
{"j":4,"side":10,"tested":184655,"admissible":117,"negative":null,"admissibleZero":null}
```

Here “admissible” means that the row/column sets admit a perfect matching in the positive-support graph (10). Every such minor was strictly positive, not merely nonzero. The check includes every minor of sizes $2$ through $10$ in the leading $10\times10$ matrix; $1\times1$ nonnegativity is immediate from (10).

This is finite evidence for TN and support-strict TN. It does not prove any all-index statement, and it does not cover $j\ge5$.

### 10.3 Concrete proof routes and the exact gaps

**Route A — support-strict TN via a planar network or bidiagonal factorization.**

Construct, for every finite truncation of $(a_{n,m})$, a planar network or a product of nonnegative bidiagonal matrices. Uniqueness of the nonintersecting-path matching would upgrade TN to support-strict TN. This would supply the missing all-minor statement directly.

*Gap:* no such network/factorization is known for the target. The natural left-production factorization is not TN, as (15)–(16) show. A successful construction must use the label-choice weights $\binom nr$ in a genuinely new way.

**Route B — direct determinant positivity from (2).**

Use
\[
a_{n,m}=n![z^n]\frac{C_j(z)^m}{(1-A_j(z))^{m+1}}
\]
and Cauchy–Binet to express a general minor as a sum of weights of nonintersecting combinatorial objects.

*Gap:* no positive Cauchy–Binet expansion is known. Positivity of individual coefficients is immediate, but positivity of arbitrary minors is exactly the unresolved statement.

**Route C — left production plus a different sufficient condition.**

Chen–Fu–Ruan Theorem 1.3 would give real-rootedness if a *different* left production matrix could be proved TN. The natural one fails.

*Gap:* find an alternative factorization whose production matrix is TN, or extend the theorem to a class that includes this array despite the negative natural-production minor.

**Route D — rootwise induction instead of full TN.**

Use (4), or a derivative-corrected version, and prove a rootwise inequality at each negative root. The issue #108 and RealRooted research notes identify this as a possible route.

*Gap:* the signed $(t-1)\binom n{j-1}H_{n-j+1}$ term needs a quantitative domination lemma. Finite Newton/Sturm tests and the two edge asymptotics do not provide that lemma for the middle roots.

Even if (11) is proved, **interlacing still needs an additional step**. A clean sufficient target would be support-strict $\mathrm{TN}_\infty$ together with row Pólya frequency and the relevant no-common-root/degree conditions, followed by the classical strict-TP/variation-diminution or Wronskian criterion. Chen–Fu–Ruan Theorem 1.3 alone gives real-rootedness, not interlacing.

---

## 11. Source list

1. **Chen–Zhu (2026), exact target enumeration.**  
   Yu Chen and Wenxiang Zhu, *A Snail Race Problem*.  
   https://arxiv.org/abs/2607.02085  
   https://arxiv.org/html/2607.02085

2. **OEIS column for $j=2,m=1$.**  
   https://oeis.org/A397265

3. **Alexandersson's current conjecture statement.**  
   https://www.symmetricfunctions.com/realRootedOpenProblems.htm#set-partitions-and-trees

4. **Alexandersson's ordered-lift theorem and proof sketch.**  
   https://www.symmetricfunctions.com/realRootedSetPartitions.htm#bigBlockSetPartitions

5. **Issue #108 and cited research commits.**  
   https://github.com/1337isnot1337/solved-conjectures/issues/108  
   https://github.com/1337isnot1337/solved-conjectures/commit/04de7069df504ff494d81085f3bbb6b0d62bb3ed  
   https://github.com/1337isnot1337/solved-conjectures/commit/e76826ae91fb6415083377ee25c5505e6462b272  
   https://github.com/1337isnot1337/solved-conjectures/blob/f1a5ac298febcfd20b19e1ef3465325eb62e18d6/research/alexandersson-ordered-set-partitions/edge.md  
   https://github.com/1337isnot1337/solved-conjectures/blob/e39bd018c0f980ff394daa7764272226890c6d14/research/alexandersson-ordered-set-partitions/left-edge.md

6. **Current RealRooted tracking issue for ordered block weights.**  
   https://github.com/PerAlexandersson/RealRooted/issues/2590

7. **Ordinary Fubini roots and current zero asymptotics.**  
   Jalowy, Kabluchko, Marynych, *Zeros and Exponential Profiles of Polynomials II: Examples*.  
   https://arxiv.org/abs/2509.11248  
   https://link.springer.com/article/10.1007/s00365-026-09775-2

8. **Proper Riordan total positivity.**  
   Xi Chen, Huyile Liang, Yi Wang, *Total positivity of Riordan arrays*.  
   https://arxiv.org/abs/1601.05637  
   https://arxiv.org/html/1601.05637

9. **Exponential Riordan total positivity.**  
   Bao-Xuan Zhu, *Total positivity from the exponential Riordan arrays*.  
   https://arxiv.org/abs/2006.14485  
   https://arxiv.org/html/2006.14485v3

10. **Left production matrices and row real-rootedness.**  
    Xi Chen, Lang Fu, Jiajie Ruan, *A unified planar network approach to total positivity of combinatorial matrices and real-rootedness of polynomials*.  
    https://arxiv.org/abs/2512.08369  
    https://arxiv.org/html/2512.08369

11. **Reciprocal-series and PF-sequence interlacing theorems.**  
    Petter Brändén and Leonardo Saud Maia Leite, *Totally nonnegative matrices, chain enumeration and zeros of polynomials*.  
    https://arxiv.org/abs/2412.06595  
    https://arxiv.org/html/2412.06595v3

12. **Finite PF-sequence criterion.**  
    Aissen, Edrei, Schoenberg, Whitney, *On the Generating Functions of Totally Positive Sequences*.  
    https://doi.org/10.1073/pnas.37.5.303

13. **Continuous Pólya-frequency kernels and transforms.**  
    Belton, Guillot, Khare, Putinar, *Totally positive kernels, Polya frequency functions, and their transforms*.  
    https://arxiv.org/abs/2006.16213  
    https://doi.org/10.1007/s11854-022-0259-7

14. **Interlacing, weighted renewal, and matrix criteria.**  
    https://www.symmetricfunctions.com/realRootedInterlacing.htm
