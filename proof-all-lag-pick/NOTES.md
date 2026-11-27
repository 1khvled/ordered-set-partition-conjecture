# All-lag Pick conjecture for the ordered big-block family

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-all-lag-pick` was created or edited.  
Status: **no proof and no counterexample found; exact structural obstruction recorded.**

## 0. Problem

Fix $j\ge2$. Let $O_n(t)=O_{n,j}(t)$ be the ordered big-block row polynomials,
with $O_0=1$, $\deg_t O_n=\lfloor n/j\rfloor$, and positive coefficients
(§1). The **all-lag Pick conjecture** is

> **(P-all)** For every $n\ge j$ and every $1\le d\le j$,
> $W(O_n,O_{n-d})(t):=O_n'(t)O_{n-d}(t)-O_n(t)O_{n-d}'(t)$ has no real zero,
> and the residues $O_n(\mu)/O_{n-d}'(\mu)$ are negative at every real root
> $\mu$ of $O_{n-d}$; equivalently $R_{n,d}:=O_n/O_{n-d}$ is a Pick/Nevanlinna
> function.

This implies the original consecutive real-rootedness/interlacing conjecture by
the IVT argument on the $d=1$ lag (Proposition B below). The $d=j$ case is the
main WRS/Wronskian target of `proof-global-wronskian/NOTES.md`.

## 1. Setup

With $A=\sum_{r=1}^{j-1}z^r/r!$, $E=\sum_{r=0}^{j-1}z^r/r!$, $C=e^z-E$,
$B=1-A$, $D=B-tC$, $F=1/D$, $O_n=n![z^n]F$, the first-block recurrences are

$$
O_n=\sum_{r=1}^{j-1}\binom nr O_{n-r}+H_n,\qquad
H_n=t\sum_{r=j}^{n}\binom nr O_{n-r},\qquad n\ge j. \tag{1}
$$

The differential form is $O_n=A_n+\binom{n-1}{j-1}tQ(O_{n-j})$, $Q(q)=q+(t-1)q'$,
$A_n=nO_{n-1}+2TO_{n-1}+\sum_{r=2}^{j-1}\binom nr O_{n-r}$ (proved in
`proof-final-wrs-corrected/NOTES.md`).

Degree: $m_n:=\deg O_n=\lfloor n/j\rfloor$. For $1\le d\le j$, writing
$n=qj+s$, $0\le s<j$,

$$
 m_n-m_{n-d}=\begin{cases}
 1,& d>s,\\
 0,& d\le s.
 \end{cases} \tag{2}
$$

Thus the Pick ratio has degree difference $1$ or $0$; both are allowed in
(P-all).

## 2. Pick equivalences (per lag)

Fix $n,d$ with $1\le d\le j$, $p=O_n$, $q=O_{n-d}$, $e=m_{n-d}$, and assume
the standing induction hypothesis that $q$ has $e$ simple negative roots
$\mu_1<\cdots<\mu_e<0$ and $p$ has positive leading coefficient. Let
$R=p/q=a t^\delta+b+\sum_{i=1}^{e}r_i/(t-\mu_i)$ with $\delta=m_n-m_{n-d}\in\{0,1\}$,
$a=\operatorname{lc}(p)/\operatorname{lc}(q)>0$, $r_i=p(\mu_i)/q'(\mu_i)$.
Then

$$
W(p,q)=q^2R'(t)=q^2\Bigl(a\delta-\sum_{i=1}^{e}\frac{r_i}{(t-\mu_i)^2}\Bigr),
\qquad W(p,q)(\mu_i)=-r_i q'(\mu_i)^2. \tag{3}
$$

**Proposition A (per-lag Pick equivalence).** Under the standing hypothesis,
the following are equivalent: (A1) $W(p,q)>0$ on $\mathbb R$; (A2) $W(p,q)>0$
on $(-\infty,0]$; (A3) $r_i<0$ for all $i$; (A4) $R=p/q$ is Pick
($\Im R(z)>0$ for $\Im z>0$); (A5) $p(\mu_i)/q'(\mu_i)<0$.

*Proof.* Identical to `proof-global-wronskian/NOTES.md`, Proposition A: (A1)
$\Rightarrow$(A2) trivial; (A2)$\Rightarrow$(A3) by (3) at $\mu_i<0$;
(A3)$\Rightarrow$(A1) because $a\delta-\sum r_i/(t-\mu_i)^2\ge a\delta\ge0$ and
is $>0$ unless $\delta=0$ and all $r_i=0$, impossible; (A3)$\Leftrightarrow$
(A4) by the Nevanlinna representation $R(z)=a\delta z+b+\sum w_i/(\mu_i-z)$,
$w_i=-r_i$; (A3)$\Leftrightarrow$(A5) by definition. $\square$

**Remark.** The $\delta=0$ case ($\deg p=\deg q$) is allowed and occurs
whenever $d\le n\bmod j$; it gives a degree-0 Pick ratio with $a=0$ in the
linear term.

## 3. All-lag Pick implies the original conjecture

**Proposition B.** Assume (P-all) for a fixed $j$ and all $n'\le n$. Then
$O_0,\ldots,O_n$ are real-rooted with simple negative roots and $O_{m-1}$
interlaces $O_m$ for all $1\le m\le n$.

*Proof.* Induct on $m$. The base cases $m<j$ are positive constants and $m=j$
is linear. Assume $O_{m-1}$ has $e$ simple negative roots. Apply Proposition A
to $(p,q)=(O_m,O_{m-1})$ ($d=1$): (P-all) gives (A3), hence $\operatorname{sgn}
O_m(\mu_i)=-\operatorname{sgn}O_{m-1}'(\mu_i)$ alternates over the roots of
$O_{m-1}$, and endpoint signs at $-\infty$ and $0$ give $e+1=\deg O_m$ roots by
IVT. $\square$

The same IVT argument also gives $O_{n-d}\interl O_n$ for every lag $d$ in
(P-all) whose degree difference is $1$; for degree difference $0$ it gives
alternation of $O_n$ at roots of $O_{n-d}$ (a Pick/Markov alternation), which
is (A3).

## 4. The all-lag Wronskian decomposition

Let $q_d:=O_{n-d}$, $1\le d\le j$. Bilinearity of $W(\cdot,q_d)$ applied to
(1), using $W(q_d,q_d)=0$, gives

$$
\boxed{\;
W(O_n,O_{n-d})=W(H_n,O_{n-d})
 +\sum_{r=1}^{d-1}\binom nr W(O_{n-r},O_{n-d})
 -\sum_{r=d+1}^{j-1}\binom nr W(O_{n-d},O_{n-r}).
\;} \tag{4}
$$

The three groups are:

- **boundary**: $W(H_n,O_{n-d})$ (cross-family, same $n$);
- **older-positive** (if $r<d$): $W(O_{n-r},O_{n-d})=V_{n-r,d-r}$, a lag
  $(d-r)$ Wronskian at index $m=n-r<n$;
- **newer-negative** (if $r>d$): $-W(O_{n-d},O_{n-r})=-V_{n-d,r-d}$, minus a
  lag $(r-d)$ Wronskian at index $m=n-d<n$.

Identity (4) is verified exactly by `verify_all_lag.py` for $j=2,\ldots,6$,
$n\le 4j$, all $d$ (0 failures).

## 5. Joint-induction obstruction (exact)

Assume the **all-lag induction hypothesis**: (P-all) holds for every $m<n$.
Then every $V_{m,k}$ appearing in (4) with $m<n$ is positive on $\mathbb R$.
The decomposition becomes

$$
W(O_n,O_{n-d})=\underbrace{W(H_n,O_{n-d})}_{\text{boundary}}
+\underbrace{\sum_{r=1}^{d-1}\binom nr V_{n-r,d-r}}_{\ge0\ \text{by IH}}
-\underbrace{\sum_{r=d+1}^{j-1}\binom nr V_{n-d,r-d}}_{\ge0\ \text{by IH, subtracted}}.
\tag{5}
$$

This has two structural consequences.

**5.1. The $d=j$ case closes only at the boundary.** For $d=j$ the negative sum
is empty, so (P-all) for lag $j$ reduces, under the all-lag IH, to
$W(H_n,O_{n-j})>0$ on $\mathbb R$. This is exactly the cross-family boundary
identified in `proof-global-wronskian/NOTES.md`, §5 and
`proof-final-wrs-corrected/NOTES.md`, §4.

**5.2. The $d<j$ cases have unavoidable negative terms.** For $d<j$, the
subtracted terms $V_{n-d,r-d}$ are themselves all-lag Pick Wronskians at lower
index; under IH they are strictly positive, so (5) contains genuine negative
contributions. Therefore the joint induction is **not closed** by sign alone:
it requires a quantitative domination inequality

$$
W(H_n,O_{n-d})+\sum_{r=1}^{d-1}\binom nr V_{n-r,d-r}
\;>\;\sum_{r=d+1}^{j-1}\binom nr V_{n-d,r-d},
\qquad d<j, \tag{6}
$$

in addition to the $d=j$ boundary condition. The all-lag hypothesis strengthens
the induction *input* but also makes the $d<j$ negative terms certifiably
present; it does not remove the tail-domination gap.

**5.3. Evidence on the boundary.** Exact Sturm counts show
$W(H_n,O_{n-d})$ has no real root for $j=2,\ldots,7$, $n\le\min(5j,32)$,
$1\le d\le j$ (562 pairs, 0 failures; `verify_all_lag.py`). Thus the boundary
is empirically Pick for all lags, but this is not proved.

## 6. Finite verification

`verify_all_lag.py` (SymPy exact arithmetic, Sturm/count_roots over $\mathbb Q$):

- all-lag Pick $W(O_n,O_{n-d})$, $j=2..9$, $n\le\min(5j,34)$, $1\le d\le j$:
  **1026 tested, 0 real roots**.
- cross-family all-lag $W(H_n,O_{n-d})$, $j=2..8$, $n\le\min(5j,30)$: **732
  tested, 0 real roots**.
- decomposition (4): **0 failures** over $j=2..6$, $n\le4j$.
- degree-difference-0 cases observed (105 in that range), confirming (2).
- 140 cases with $d<j$ having at least one subtracted IH-positive term in (5).
- all tested $W(O_n,O_{n-d})$ have positive leading coefficient.

These are exact finite checks, not a proof.

## 7. Precise gap

To prove (P-all) by induction on $n$ it is enough to prove, under the all-lag IH
for all $m<n$, the family of inequalities

$$
\boxed{\;
\begin{aligned}
& W(H_n,O_{n-j})>0, && d=j,\\[2pt]
& W(H_n,O_{n-d})+\!\sum_{r<d}\binom nr V_{n-r,d-r}
  > \!\sum_{r>d}\binom nr V_{n-d,r-d}, && 1\le d<j,
\end{aligned}
\;} \tag{7}
$$

for all real $t$. Equivalently, by Proposition A, prove the corresponding
residue/Bezoutian domination at every root of $O_{n-d}$.

What is **not** obtained: (7). The all-lag framing exposes a real obstruction
rather than a proof: the joint induction is sign-closed only at the main lag
$d=j$ and only after the cross-family boundary is supplied; for $d<j$ it is
quantitatively obstructed by IH-positive subtracted terms. No termwise-positive
kernel, no common-interlacer certificate, and no ratio-monotonicity identity
closing (7) was found. The $d=j$ boundary is the same WRS gap as in
`proof-global-wronskian/NOTES.md`.

## 8. Files

- `NOTES.md` — this report.
- `verify_all_lag.py` — exact decomposition and Sturm/root-count checks.

Run: `py -3.14 verify_all_lag.py` from this directory. No packages installed,
no builds, no large caches; total new data well under 1 MB.