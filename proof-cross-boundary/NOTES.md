# Cross-family boundary Wronskian attack

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-cross-boundary` was created or edited.  
Status: **no unconditional all-index proof of $W(H_{n,j},O_{n-d,j})>0$ on $\mathbb R$ was obtained.** A clean reduction is proved, the $d=j$ boundary is isolated as an exact anti-Pick/domination inequality, and the termwise/inner-family routes are shown to be obstructed. No counterexample was found in exact finite checks.

## 0. Target

Fix $j\ge2$. Let $O_n(t)=O_{n,j}(t)$ be the ordered big-block row polynomials,
$A=\sum_{r=1}^{j-1}z^r/r!$, $C=e^z-\sum_{r=0}^{j-1}z^r/r!$, $B=1-A$,
$D=B-tC$, $F=1/D$, $O_n=n![z^n]F$. Let
$$
 g_j=\frac{C}{1-A},\qquad H(z,t)=\frac{1}{1-tg_j(z)}=\sum_{n\ge0}H_n(t)\frac{z^n}{n!},
 \qquad H_n=tS_n,\qquad S_n:=\sum_{k=0}^{n-j}\binom{n}{k}O_k.
\tag{0.1}
$$
Equivalently $H_n=t\sum_{r=j}^{n}\binom nr O_{n-r}$ and
$O_n=\sum_{r=1}^{j-1}\binom nr O_{n-r}+H_n$ for $n\ge j$. The target is
$$
\boxed{\;W(H_n,O_{n-d})(t)=H_n'(t)O_{n-d}(t)-H_n(t)O_{n-d}'(t)>0
\quad (1\le d\le j,\ n\ge j,\ t\in\mathbb R),\;}
\tag{0.2}
$$
with special emphasis on $d=j$. (All derivatives are in $t$.)

## 1. Basic facts used

Degree and roots. $\deg_t O_n=\lfloor n/j\rfloor=:m_n$, $\deg_t H_n=m_n$ for
$n\ge j$, and both have positive leading coefficients. For $q_d:=O_{n-d}$,
$\deg q_d=m_{n-d}\in\{m_n-1,m_n\}$ and
$$
 m_n-m_{n-d}=\begin{cases}1,&d>s,\\0,&d\le s,\end{cases}\qquad n=qj+s,\ 0\le s<j.
\tag{1.1}
$$

Standing induction hypothesis (IH). For the residue/Pick reduction below we
assume, for all $m<n$, that $O_m$ has $m_m=\lfloor m/j\rfloor$ simple negative
roots and the all-lag Pick condition $W(O_m,O_{m-k})>0$ on $\mathbb R$ for
$1\le k\le j$, $m-k\ge0$. By Proposition A of
`proof-global-wronskian/NOTES.md` and `proof-all-lag-pick/NOTES.md`, this is
equivalent to negative residues $O_m(\mu)/O_{m-k}'(\mu)<0$ at roots $\mu$ of
$O_{m-k}$.

The product rule (problem hint) is
$$
 W(tp,q)=p q + tW(p,q).
\tag{1.2}
$$

The chapter decomposition (problem hint) is $H=\sum_{m\ge1}t^m g_j^m$ and
$F=\frac{1}{1-A_j}H$, hence
$$
 O_m=\sum_{k=0}^{m}\binom{m}{k}d_{m-k}H_k,\qquad d_r:=r![z^r]\frac{1}{1-A_j}\ge0.
\tag{1.3}
$$

## 2. Product-rule / residue reduction

Let $q_d=O_{n-d}$, $1\le d\le j$, $n\ge j$. Since $H_n=tS_n$,
$$
 W(H_n,q_d)=S_n q_d + t W(S_n,q_d).
\tag{2.1}
$$
At a root $\mu$ of $q_d$, $q_d(\mu)=0$ and therefore
$$
 W(H_n,q_d)(\mu)=\mu\, W(S_n,q_d)(\mu)=-\mu\, S_n(\mu)q_d'(\mu).
\tag{2.2}
$$
Under IH, $\mu<0$ and $q_d'(\mu)\neq0$, so
$$
 W(H_n,q_d)(\mu)>0 \iff W(S_n,q_d)(\mu)<0 \iff \frac{S_n(\mu)}{q_d'(\mu)}>0
 \iff \frac{H_n(\mu)}{q_d'(\mu)}<0.
\tag{2.3}
$$

**Proposition 2.1 (Pick equivalence for the boundary).** Assume IH for $q_d$:
$q_d$ has simple negative roots and $\deg H_n-\deg q_d\in\{0,1\}$, with
positive leading coefficients. Then the following are equivalent:
- (B1) $W(H_n,q_d)>0$ on $\mathbb R$;
- (B2) $W(H_n,q_d)>0$ on $(-\infty,0]$;
- (B3) $H_n(\mu)/q_d'(\mu)<0$ at every real root $\mu$ of $q_d$;
- (B4) $S_n(\mu)/q_d'(\mu)>0$ at every real root $\mu$ of $q_d$;
- (B5) $H_n/q_d$ is a Pick/Nevanlinna function (degree-diff 1 case), or an
  anti-Pick/monotone rational function with negative residues in the
  degree-diff 0 case in the same $W=q^2R'$ sense.

*Proof.* Exactly Proposition A of `proof-global-wronskian/NOTES.md`, applied to
$p=H_n,q=q_d$, $\delta=\deg H_n-\deg q_d\in\{0,1\}$: $W=q^2(a\delta-\sum r_i/(t-\mu_i)^2)$, $r_i=H_n(\mu_i)/q_d'(\mu_i)$, and $W(\mu_i)=-r_i q_d'(\mu_i)^2$. The equivalence (B3)$\Leftrightarrow$(B4) is (2.2) with $\mu_i<0$. $\square$

Thus the target (0.2) is an exact residue/Bezoutian domination problem, as
requested.

## 3. The $d=j$ boundary is clean

Let $q:=q_j=O_{n-j}$, $e=m_n-1=\deg q$, $m=n-j$. The sum $S_n$ has the same top
term as $\binom{n}{j}q$:
$$
 S_n=\binom{n}{j}q + R_n,\qquad R_n:=\sum_{k=0}^{m-1}\binom{n}{k}O_k.
\tag{3.1}
$$
Since $W(q,q)=0$,
$$
 W(S_n,q)=\sum_{k=0}^{m-1}\binom{n}{k}W(O_k,q)
         = -\sum_{k=0}^{m-1}\binom{n}{k}W(q,O_k).
\tag{3.2}
$$

**Proposition 3.1 ($d=j$ reduction).** Under IH for $q=O_{n-j}$,
$$
 W(H_n,q)>0 \ \text{on }\mathbb R
 \quad\Longleftrightarrow\quad
 W(S_n,q)<0 \ \text{on }\mathbb R
 \quad\Longleftrightarrow\quad
 \sum_{k=0}^{n-j-1}\binom{n}{k}W(O_{n-j},O_k)>0\ \text{on }\mathbb R.
\tag{3.3}
$$

*Proof.* $S_n$ and $q$ have equal degree $e$ and positive leading coefficients,
so $R=S_n/q$ has degree difference 0 and $W(S_n,q)=q^2 R'$. No real roots and
negative value (equivalently $R'<0$) is equivalent to all residues
$S_n(\mu)/q'(\mu)>0$. By (2.3) this is equivalent to $W(H_n,q)(\mu)>0$ at
roots of $q$, and Proposition 2.1 upgrades it to all $\mathbb R$. The last
equivalence is (3.2). $\square$

**Consequence (stronger than the main WRS at residues).** Inserting
$O_n=H_n+\sum_{r=1}^{j-1}\binom nr O_{n-r}$ into $W(O_n,q)$ gives
$$
 W(O_n,q)=W(H_n,q)+\sum_{r=1}^{j-1}\binom nr W(O_{n-r},q).
\tag{3.4}
$$
For $r<j=d$, IH gives $W(O_{n-r},O_{n-j})>0$. Hence $W(H_n,q)>0$ implies
$W(O_n,q)>0$, but not conversely: at a root $\mu$ of $q$,
$$
 \frac{H_n(\mu)}{q'(\mu)}=\frac{O_n(\mu)}{q'(\mu)}+\sum_{r=1}^{j-1}\binom nr
 \left| \frac{O_{n-r}(\mu)}{q'(\mu)}\right|.
\tag{3.5}
$$
Thus the $d=j$ boundary (0.2) is a stronger residue inequality than the main
$W(O_n,O_{n-j})>0$ boundary.

## 4. Why the $d=j$ domination is nontrivial

Write $m=n-j$ and split (3.3) into close and far lags:
$$
 \Sigma_n:=\sum_{k=0}^{m-1}\binom{n}{k}W(O_m,O_k)=
 \underbrace{\sum_{k=\max(0,m-j)}^{m-1}\binom{n}{k}W(O_m,O_k)}_{\Sigma_n^{\mathrm{close}}}
 +\underbrace{\sum_{k=0}^{m-j-1}\binom{n}{k}W(O_m,O_k)}_{\Sigma_n^{\mathrm{far}}}.
\tag{4.1}
$$

For close lags $1\le m-k\le j$, IH gives $W(O_m,O_k)>0$ on $\mathbb R$; hence
$\Sigma_n^{\mathrm{close}}>0$ termwise. For far lags $m-k>j$, degree difference
$\deg O_m-\deg O_k$ can be $\ge2$, so $W(O_m,O_k)$ cannot be globally
Pick/positive on $\mathbb R$ and indeed has real roots in tests
(`probe_lag.py`: e.g. $j=3,n=12,m=9$, lags 4--9 all have real-root Wronskians).
Moreover the residue of the reverse ratio at a root $\mu$ of $O_m$,
$W(O_m,O_k)(\mu)=O_m'(\mu)O_k(\mu)$, is not sign-definite for far lags
(`probe_root_sign.py`: many negative values for lags $>j$). Therefore the
positivity of $\Sigma_n$ is a genuine domination of IH-positive close-lag
terms over indefinite far-lag terms, not a termwise consequence of the all-lag
IH.

## 5. The $d<j$ boundary

For $1\le d<j$, use the full decomposition (proof-all-lag-pick (4)) at a root
$\mu$ of $q_d=O_{n-d}$:
$$
 \frac{H_n(\mu)}{q_d'(\mu)}=
 \frac{O_n(\mu)}{q_d'(\mu)}
 -\sum_{r=1}^{d-1}\binom nr\frac{O_{n-r}(\mu)}{q_d'(\mu)}
 -\sum_{r=d+1}^{j-1}\binom nr\frac{O_{n-r}(\mu)}{q_d'(\mu)}.
\tag{5.1}
$$
IH gives $O_{n-r}(\mu)/q_d'(\mu)<0$ for $r<d$ and $>0$ for $r>d$. Hence
$$
 W(H_n,q_d)(\mu)>0 \iff \frac{O_n(\mu)}{q_d'(\mu)} < 
 \sum_{r=1}^{d-1}\binom nr\frac{O_{n-r}(\mu)}{q_d'(\mu)}
 +\sum_{r=d+1}^{j-1}\binom nr\frac{O_{n-r}(\mu)}{q_d'(\mu)}.
\tag{5.2}
$$
On the right, the $r<d$ IH-terms are negative and the $r>d$ IH-terms are
positive, so the bound is a mixed-sign quantitative condition. In particular,
for $d<j$ the boundary $W(H_n,q_d)>0$ is neither implied by nor does it imply
$W(O_n,q_d)>0$ alone without that quantitative bound. This is the same
obstruction as `proof-all-lag-pick/NOTES.md` ��5.2, now seen at residue level.

## 6. Chapter-decomposition identity and inner-family obstruction

Using (1.3) in $W(H_n,O_m)$ gives the chapter decomposition
$$
 W(H_n,O_m)=\sum_{k=0}^{m}\binom{m}{k}d_{m-k}W(H_n,H_k),\qquad d_r\ge0.
\tag{6.1}
$$
This is a positive-coefficient combination of inner Wronskians. If the inner
family had all-lag Pick, $W(H_n,H_k)>0$, (6.1) would imply (0.2). However,
exact Sturm counts show this is false: $W(H_n,H_{n-j})$ has real roots for all
tested $j=2,\ldots,6$, $n\le\min(5j,30)$ (65 tested, 65 with real roots;
`probe_HH.py`), and most $W(H_n,H_k)$ have real roots (`probe_more.py`).
Thus (6.1) is a true identity but not a termwise-positive proof route.

## 7. Double-resolvent kernel

Let $S_n=n![z^n](C(z)F(z))$, $O_m=m![w^m]F(w)$, $\Phi=C(z)B(w)-C(w)B(z)$,
$\Psi=D(z)D(w)+t\Phi(z,w)=B(z)D(w)-tC(w)D(z)$. The same Cauchy extraction as
`proof-global-wronskian/NOTES.md` Proposition D gives
$$
 W(S_n,O_{n-d})=n!(n-d)![z^n w^{n-d}]
 \frac{C(z)\Phi(z,w)}{D(z)^2 D(w)^2},
\tag{7.1}
$$
$$
 W(H_n,O_{n-d})=n!(n-d)![z^n w^{n-d}]
 \frac{C(z)\Psi(z,w)}{D(z)^2 D(w)^2}.
\tag{7.2}
$$
Identity (7.2) is verified symbolically for $j=2,n=6,d=2$ in
`check_kernel2.py`. The $C(z)$ factor distinguishes the boundary kernel from
the main $W(O_n,O_{n-j})$ kernel $\Phi/(D^2D^2)$; $\Phi$ remains indefinite on
$(0,\infty)^2$ (Prop. D, loc. cit.), so no termwise-positive real-contour
kernel proof is available from this representation.

## 8. Finite verification

All checks use exact SymPy/QQ arithmetic.

- $W(H_n,O_{n-d})$ real-root count: $j=2..7$, $n\le\min(5j,30)$, $1\le d\le j$:
  548 tested, 0 with real roots (`probe.py`).
- $S_n(\mu)/O_{n-d}'(\mu)>0$ at roots of $O_{n-d}$: $j=2..6$, $n\le\min(4j,24)$:
  435 residue tests, 0 failures (`probe_S.py`).
- $W(H_n,O_{n-d})$ coefficientwise nonnegative: $j=2..7$, $n\le\min(5j,28)$:
  522 tested, 0 failures (`probe_coeff.py`).
- $W(S_n,O_{n-j})$ coefficientwise nonpositive and real-root free:
  $j=2..7$, $n\le\min(5j,28)$: 78 tested, 0 failures (`probe_more.py`,
  `probe_WS.py`, `probe_WS_sign.py`).
- Far-lag $W(O_{n-j},O_k)$ real roots: e.g. $j=3,n=12$, lags 4--9 have real
  roots (`probe_lag.py`).
- Reverse-residue sign failures $W(O_{n-j},O_k)(\mu)=O_{n-j}'(\mu)O_k(\mu)\le0$
  at roots $\mu$ of $O_{n-j}$ for far lags: 1625 tested, 404 failures
  (`probe_root_sign.py`; expected for large lags).
- $W(H_n,H_{n-j})$ real-root obstruction: $j=2..6$, $n\le\min(5j,30)$: 65
  tested, 65 with real roots (`probe_HH.py`).
- Chapter identity (1.3) checked exactly for $j=3$, $n\le15$ (`check_dH.py`).

## 9. Exact obstruction / remaining gap

The $d=j$ boundary is exactly equivalent, under IH, to the anti-Pick/domination
inequality
$$
\boxed{\;\sum_{k=0}^{n-j-1}\binom{n}{k}W(O_{n-j},O_k)(t)>0,
\qquad \forall t\in\mathbb R,\;}
\tag{9.1}
$$
or, at roots $\mu$ of $O_{n-j}$, to $\sum_{k=0}^{n-j-1}\binom{n}{k}
O_k(\mu)/O_{n-j}'(\mu)>0$. The close-lag part is IH-positive termwise; the
far-lag part is not (real roots and reverse-residue sign failures occur), so
termwise closure is false. The chapter decomposition (6.1) gives a positive
combination of inner Wronskians, but those inner Wronskians themselves fail
Pick. The kernel (7.2) is indefinite in the sense of Prop. D.

For $d<j$, the boundary is the stronger mixed-sign residue domination (5.2).
No common-interlacer, no termwise-positive kernel, and no Bezoutian domination
certificate closing (9.1)/(5.2) was found. A proof of (0.2) would imply the
main $d=j$ WRS/Wronskian boundary by (3.4) and so is at least as hard as the
existing global Wronskian gap; this work does not close it.
