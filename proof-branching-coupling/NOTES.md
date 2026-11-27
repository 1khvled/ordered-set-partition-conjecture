# Measure-valued branching coupling attempt for the cross-boundary Pick inequality

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-branching-coupling` was created or edited.  
Outcome: **no proof of $W(H_{n,j},O_{n-d,j})>0$ was obtained.** A natural
measure-valued branching / spine coupling is constructed (Theorems 1--3): it
realizes the chapter renewal as a positive CB-type Möbius flow with an explicit
size-biased spine martingale, and it identifies the boundary residue sum in (G)
as a **signed cross-variation / renewal-mass functional**. The construction makes
the obstruction precise in two stronger forms: the coefficient-extraction
functional is non-positive on the Pick cone of the flow (already C1 of
`proof-chapter-renewal/NOTES.md`), and the signed cross-variation cannot be made
nonnegative by any coupling that respects the chapter spine because the
far-lag boundary mass has no global sign and the reciprocal kernel $1/g=B/C$ is
not a Stieltjes/Pick kernel on $(0,\infty)$ (Theorems 4--5). The exact remaining
gap is recorded as a signed domination inequality (G').

---

## 0. Setup and target

Fix $j\ge2$ and use the chapter-renewal dictionary:
$$
A=\sum_{r=1}^{j-1}\frac{z^r}{r!},\quad C=e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!},
\quad B=1-A,\quad p=\frac{z^{j-1}}{(j-1)!},\quad g=\frac{C}{B},
$$
$$
D=B-tC,\quad F=1/D=\sum_{n\ge0}O_n(t)\frac{z^n}{n!},\quad
H=B F=1/(1-tg)=\sum_{n\ge0}H_n(t)\frac{z^n}{n!}.
$$
Identities (all proved in cited notes): $A'=1+A-p$, $C'=C+p$, $B'=B-2+p$,
$B g'=(2-p)g+p$, $F'=(2+(t-1)p)F^2-F$.

The boundary target is, for $1\le d\le j$, $q:=O_{n-d}$,
$$
W(H_n,q)(t)=H_n'(t)q(t)-H_n(t)q'(t)>0,\qquad t\in\mathbb R. \tag{0.1}
$$
Under the standing simple-negative-root/all-lag Pick IH for $q$,
Proposition 2.1 of `proof-cross-boundary/NOTES.md` and Theorem B2 of
`proof-chapter-renewal/NOTES.md` give the equivalent residue form
$$
\boxed{\;W(H_n,q)>0\text{ on }\mathbb R \iff G_q(\mu):=\sum_{r=j}^{n}\binom nr
\frac{O_{n-r}(\mu)}{q'(\mu)}>0\ \ \text{for every root }\mu<0\text{ of }q.\;}
\tag{G}
$$
Section 5 of `proof-chapter-renewal/NOTES.md` asks for a measure-valued
martingale or branching coupling whose cross-variation at the poles of $q$ is
exactly the left side of (G), with strictly positive value.

---

## 1. Chapter branching measure

### 1.1 The renewal/branching random object

Let $\mathfrak C_n$ be the set of ordered chapter compositions of $[n]$:
an element is a sequence $\pi=(C_1,\ldots,C_m)$ of disjoint nonempty blocks whose
union is $[n]$, in order, where every chapter $C_i$ has the form
"a marked block of size $\ge j$ followed by a string of unmarked blocks of
sizes in $\{1,\ldots,j-1\}$". Equivalently, by the renewal equation
$H=1+tgH$, the chapter EGF is $g(z)=C/B=\sum_{k\ge j}a_k z^k/k!$ with
$a_k=k![z^k]g\ge0$ and $a_k=0$ for $k<j$.

Define the positive finite measure on $\mathfrak C_n$ by
$$
\mathsf M_{n,t}(\{\pi\}):=t^{m(\pi)}\prod_{i=1}^{m(\pi)} a_{|C_i|},
\qquad m(\pi)=\text{number of chapters}. \tag{1.1}
$$
(Combinatorially this is the ordinary weight $t^m\prod a_{|C_i|}/|C_i|!$
times the labeling factor $n!/\prod |C_i|!$; we keep the $n!$-normalized
convention matching $H_n=n![z^n]H$.) Then
$$
H_n(t)=\sum_{\pi\in\mathfrak C_n}\mathsf M_{n,t}(\{\pi\}),\qquad
H_n(t)=\delta_{n,0}+t\sum_{r=j}^{n}\binom nr a_r H_{n-r}(t). \tag{1.2}
$$

For $t>0$, normalize $\mathbb P_{n,t}:=\mathsf M_{n,t}/H_n(t)$. Let $M=M(\pi)$ be
the chapter count and $R_1=|C_1|$ the first chapter size (when $n\ge j$). Then
$$
\mathbb P_{n,t}(R_1=r)=\frac{t\binom nr a_r H_{n-r}(t)}{H_n(t)},\qquad r=j,\ldots,n.
\tag{1.3}
$$

### 1.2 Size-biased spine martingale (positive $t$)

For $t>0$, $n\ge j$, the identity $H_n=t\sum_{r=j}^n\binom nr a_r H_{n-r}$
gives the one-step spine decomposition: conditionally on the first chapter size
$R_1=r$, the remainder is an independent $\mathbb P_{n-r,t}$ chapter composition.
The size-biased (Palm) measure by chapter count is
$$
\widehat{\mathbb P}_{n,t}(A):=\frac{\mathbb E_{n,t}[M\,1_A]}{\mathbb E_{n,t}[M]}.
\tag{1.4}
$$
Since $tH_n'(t)=\mathbb E_{n,t}[M]H_n(t)$, (1.4) is a positive martingale tilt.
The spine recursion is
$$
\mathbb E_{n,t}[M]=1+\sum_{r=j}^{n}\mathbb P_{n,t}(R_1=r)\,\mathbb E_{n-r,t}[M],
\tag{1.5}
$$
and differentiating (1.2) gives the Wronskian-as-covariance identity
$$
W(H_n,q)(t)=\frac{H_n(t)q(t)}{t}\Bigl(\mathbb E_{n,t}[M]-\mathbb E^q_{t}[M_q]\Bigr),
\qquad t>0,\ q(t)>0, \tag{1.6}
$$
where $\mathbb E^q_t$ is the analogous chapter-count mean for $q=O_{n-d}$ after
its own renewal decomposition. Identity (1.6) is exact for $t>0$ but uses
positive normalization; it does not extend to the negative roots $\mu<0$ of $q$
where the target (G) lives.

---

## 2. Continuous-state branching (CB) coupling of the Möbius flow

### 2.1 Quadratic CB mechanism

On $I_B:=(0,z_B)$, define the inhomogeneous quadratic branching mechanism
$$
\Psi_z(\lambda):=\gamma(z)\lambda+\beta(z)\lambda^2,\qquad
\gamma(z)=\frac{2-p(z)}{B(z)},\qquad \beta(z)=\frac{p(z)}{B(z)}.
\tag{2.1}
$$
Since $B>0$, $p=z^{j-1}/(j-1)!>0$ on $I_B$, we have $\beta>0$. Moreover
$2-p>0$ on $I_B$ in all tested cases $j=2,\ldots,9$ (verified in
`verify_coupling.py`, command in §8 notes), so $\gamma>0$ on the germ.  The
Riccati equation of Theorem A is exactly the CB flow equation
$$
\partial_z\phi_z(t_0)=-\Psi_z(\phi_z(t_0)),\qquad \phi_0(t_0)=t_0.
\tag{2.2}
$$

Let $(X_z)_{z\in I_B}$ be the inhomogeneous continuous-state branching (CB)
process with quadratic mechanism $\Psi_z$, a nonnegative real-valued branching
process (the one-dimensional measure-valued branching case). By the standard
quadratic CB Riccati semigroup / Lamperti-transform construction, it exists as a
nonnegative strong Markov process with Laplace functional
$$
\mathbb E_x\!\left[e^{-t_0 X_z}\right]=e^{-x\phi_z(t_0)},\qquad x\ge0,\ t_0\ge0.
\tag{2.3}
$$
Thus $\phi_z$ is the **branching coupling flow**: it couples the initial mass
$t_0$ at $z=0$ to the mass at $z$ by a positive measure-valued branching
mechanism.

### 2.2 Chapter EGF as an affine observable of the CB flow

The transport identity (A4) of `proof-chapter-renewal/NOTES.md` is
$$
H(z,\phi_z(t_0))=1+e^{-z}C(z)t_0=:A_0(z)+A_1(z)t_0,\qquad A_1(z)=e^{-z}C(z)>0.
\tag{2.4}
$$
Equivalently, evaluating the chapter EGF along the CB flow gives a positive
affine observable of the initial mass:
$$
H(z,\phi_z(t_0))=A_0(z)+A_1(z)\,t_0,\qquad A_1(z)=e^{-z}C(z)>0.
\tag{2.5}
$$
The CB process provides the positive flow $\phi_z$ transporting the mass
variable $t_0$; the chapter observable is the affine first integral (2.5). This
is the sense in which the chapter renewal is a continuous-state branching
coupling at EGF level: the $z$-evolution is the positive CB Riccati flow, and
the coefficient-level boundary target is the signed cross-variation of that flow
at poles of $O_{n-d}$.

### 2.3 Coefficient extraction is signed

The coefficient functional $\Lambda_n(f):=n![z^n]f(z)$ satisfies
$H_n(t)=\Lambda_n(H(\cdot,t))$. Cauchy's formula represents $\Lambda_n$ as a
signed/complex contour functional, not a positive measure on $I_B$. Theorem C1
of `proof-chapter-renewal/NOTES.md` shows this non-positivity is essential:
no positive measure representation of $H_n$ can hold for $n\ge2j$.  Hence a
positive CB coupling at the EGF level does not automatically pass to a positive
coefficient-level coupling.

---

## 3. Signed cross-variation coupling at boundary poles

This section constructs the object requested by §5 of
`proof-chapter-renewal/NOTES.md`: a signed martingale/coupling whose
cross-variation at a pole $\mu$ of $q=O_{n-d}$ is exactly $G_q(\mu)$.

### 3.1 Spine decomposition of the chapter renewal

For $t>0$, let $\mathbb P_{n,t}$ be the chapter measure (1.3).  Introduce the
first-chapter spine random variable $R\in\{j,\ldots,n\}$ with
$$
\mathbb P_{n,t}(R=r)=\frac{t\binom nr a_r H_{n-r}(t)}{H_n(t)},
\qquad \sum_{r=j}^n \mathbb P_{n,t}(R=r)=1. \tag{3.1}
$$
Conditionally on $R=r$, the remaining composition has law $\mathbb P_{n-r,t}$.
This is the standard renewal spine coupling: it couples $H_n$ to its
substructures $H_{n-r}$ through a positive probability kernel.

The same decomposition at the EGF level is the CB transport (2.4): the $z$-flow
couples the initial mass $t_0$ to the full chapter mass.  Coefficient
extraction turns the CB flow into the discrete spine (3.1).

### 3.2 Signed residue mass at a pole

Let $\mu<0$ be a simple root of $q=O_{n-d}$, $q'(\mu)\ne0$. Define the signed
residue mass on chapter sizes by
$$
\Xi_{n,d,\mu}(r):=\binom nr \frac{O_{n-r}(\mu)}{q'(\mu)},\qquad r=j,\ldots,n.
\tag{3.2}
$$
No positivity of $\Xi_{n,d,\mu}$ is assumed. Its total mass is
$$
G_q(\mu)=\sum_{r=j}^{n}\Xi_{n,d,\mu}(r). \tag{3.3}
$$

### 3.3 Cross-variation identity

Using $H_n=t\sum_{r=j}^{n}\binom nr O_{n-r}$ and $W(tf,q)=fq+tW(f,q)$,
$$
W(H_n,q)=\sum_{r=j}^{n}\binom nr\left(O_{n-r}q+tW(O_{n-r},q)\right). \tag{3.4}
$$
At $q(\mu)=0$, $W(O_{n-r},q)(\mu)=-O_{n-r}(\mu)q'(\mu)$, hence
$$
W(H_n,q)(\mu)=-\mu q'(\mu)\sum_{r=j}^{n}\binom nr O_{n-r}(\mu),
\qquad
\frac{W(H_n,q)(\mu)}{q'(\mu)^2}=-\mu\,G_q(\mu). \tag{3.5}
$$
Since $\mu<0$, the boundary Pick condition is
$$
\boxed{\;W(H_n,q)(\mu)>0\iff \mathcal X_{n,d}(\mu):=-\mu G_q(\mu)>0\iff G_q(\mu)>0.\;}
\tag{3.6}
$$

### 3.4 The constructed coupling/martingale

The construction requested by the problem is the following signed coupling.

- **Probability/flow layer ($z$).** The CB process (2.3) with mechanism (2.1)
  gives a positive measure-valued branching coupling for the EGF flow
  $\phi_z$.
- **Spine layer (discrete $n$).** For $t>0$, the spine kernel (3.1) is a
  positive martingale decomposition of $H_n$ into chapter submasses.
- **Boundary layer (pole $\mu$).** Evaluate the spine mass at a pole of $q$ by
  the signed residue functional $f\mapsto f(\mu)/q'(\mu)$. The resulting
  signed measure is $\Xi_{n,d,\mu}$ in (3.2). Its cross-variation with the
  Wronskian functional is exactly $-\mu G_q(\mu)$ by (3.5).

Thus the cross-variation at the pole is constructed exactly as requested.  The
missing property is not existence/identification of the cross-variation but its
positivity: $\Xi_{n,d,\mu}$ is signed, and no positive domination has been
proved.

### 3.5 Double-resolvent kernel form

Equivalently, by (7.2) of `proof-cross-boundary/NOTES.md`,
$$
W(H_n,O_{n-d})=n!(n-d)![z^n w^{n-d}]\frac{C(z)\Psi(z,w)}{D(z,t)^2D(w,t)^2},
\quad \Psi=B(z)D(w,t)-tC(w)D(z,t). \tag{3.7}
$$
A positive bivariate kernel/coupling proof would require a positive
representation of this coefficient.  The kernel is indefinite on the relevant
real domain (Proposition D of `proof-global-wronskian/NOTES.md`), so (3.7)
agrees with the signedness obstruction in §4.

---

## 4. Obstruction A: non-positivity of coefficient extraction

### Theorem 4.1 (coefficient extraction is not a positive flow functional)

There is no finite positive measure $\mu$ supported in $(0,z_B)$ and no positive
kernel $K_n(z,\cdot)$ on $(-\infty,0]$ such that
$$
H_n(t)=\int K_n(z,t)\,H(z,t)\,d\mu(z),\qquad K_n(z,t)\ge0,
\tag{4.1}
$$
and such that $K_n$ preserves the Wronskian sign needed in (0.1). In particular,
the Möbius flow positivity (2.2) cannot be transferred to (G) by positive
integration over the flow time $z$.

*Proof.* The first clause is Theorem C1 of `proof-chapter-renewal/NOTES.md`:
for $n\ge2j$, $\deg H_n\ge2$ while any positive integral of Pick slices
$H(z,\cdot)$ is a Pick function of at most linear growth. The Wronskian clause
follows because a positive kernel preserving pointwise sign cannot manufacture a
degree-2 numerator from degree-$\le1$ Pick slices; more directly, (3.5) shows the
actual coefficient kernel is the indefinite double-resolvent kernel, not a
positive one. $\square$

### Corollary 4.2

Any measure-valued branching coupling proving (G) must use either a
**signed** coefficient measure, a **signed** martingale, or a coupling whose
positivity is proved after coefficient extraction by a domination argument; it
cannot be a positive-kernel pushforward of the EGF flow alone.

---

## 5. Obstruction B: reciprocal-kernel / renewal obstruction

The second suggested equivalent form is "a renewal theorem for the reciprocal
kernel $g(z)$ that does not require $g$ to be PF/CZDS".

### 5.1 Reciprocal kernel

We have $1/g(z)=B(z)/C(z)=\alpha(z)$. The EGF Pick pole is at
$t=\alpha(z)$. On $(0,z_B)$,
$$
\alpha'(z)=\frac{B'C-BC'}{C^2}=\frac{(B-2+p)C-B(C+p)}{C^2}
= -\frac{2C+p(C-B)}{C^2}<0, \tag{5.1}
$$
since $C>0$ and $C-B=e^z-1>0$. However, globally $\alpha$ is not monotone:
$\alpha(z_B)=0$, $\alpha(z)\sim -z^{j-1}e^{-z}/(j-1)!\to0^-$, so $\alpha$ has a
minimum on $(z_B,\infty)$ (Theorem C2 of `proof-chapter-renewal/NOTES.md`).

### 5.2 Consequence for renewal/Pick

A Stieltjes/Pick reciprocal kernel on $(0,\infty)$ would require a global sign
condition (complete monotonicity / Nevanlinna-Pick kernel monotonicity) on
$\alpha$ or on the resolvent. Since $\alpha$ reverses monotonicity, no global
Pick/Stieltjes renewal theorem applies to $1/g$ on the full domain. On the small
germ $(0,z_B)$ the kernel is oriented, but Corollary 4.2 blocks passing that
germ positivity through coefficient extraction.

### Theorem 5.3 (reciprocal-kernel obstruction)

The reciprocal kernel $1/g=B/C$ is strictly decreasing on $(0,z_B)$ but not on
$(0,\infty)$; therefore it does not satisfy the global total-positivity /
Stieltjes monotonicity hypothesis required by the standard renewal
domination theorems that would imply (G). Any renewal theorem proving (G) must
operate after coefficient extraction, at the level of the signed residue sum
$G_q(\mu)$, not at the level of the global EGF kernel.

*Proof.* (5.1) gives decrease on $(0,z_B)$; the boundary/asymptotic values give
a minimum on $(z_B,\infty)$. Theorem C2 gives the resulting sign reversal of the
$2\times2$ resolvent minor. A global Stieltjes/Pick renewal domination would
force a fixed minor sign; contradiction. The final clause is Corollary 4.2.
$\square$

---

## 6. Obstruction C: the signed cross-variation is not termwise positive

Even after coefficient extraction, the cross-variation $\mathcal X_{n,d}(\mu)
=-\mu G_q(\mu)$ is a signed sum. Write $q=O_{n-d}$ and split by lag
$\ell=r-(n-d)=d - k$, $k=n-r$:
$$
G_q(\mu)=\sum_{\ell=1}^{d-1}\binom n{n-d+\ell}\frac{O_{n-d+\ell}(\mu)}{q'(\mu)}
+\sum_{\ell=d}^{n-(n-d)}\binom n{n-d+\ell}\frac{O_{n-d+\ell}(\mu)}{q'(\mu)}.
\tag{6.1}
$$
Equivalently, in the close/far lag form of `proof-cross-boundary/NOTES.md` §4,
the older terms (relative to $q$ in the all-lag decomposition) are IH-controlled
while far-lag terms have no fixed sign. In the renewal form (B2),
$O_{n-r}(\mu)/q'(\mu)= -W(O_{n-r},q)(\mu)/q'(\mu)^2$ can be positive or
negative for far lags (Theorem B witness: $j=2,n=8,q=O_6,r=5$ gives a far-lag
negative Wronskian coefficient; residue-level far-lag sign failures are recorded
in `proof-cross-boundary/NOTES.md` §8).

Thus no coupling that proves nonnegativity termwise over chapters/lags can prove
(G). A valid coupling/martingale proof must produce a **global domination**
inequality:
$$
\boxed{\;\sum_{r=j}^{n}\binom nr \frac{O_{n-r}(\mu)}{O'_{n-d}(\mu)}>0
\quad (\mu\in\mathrm{Roots}(O_{n-d})\cap(-\infty,0)).\;} \tag{G'}
$$
This is exactly the same gap as (G), now certified as irreducible to
termwise positivity, positive-kernel transfer, or global reciprocal-kernel
monotonicity.

---

## 7. What a constructed coupling would still need

The constructed objects above do give a precise target:

1. **Positive $z$-flow.** The Möbius/CB Riccati flow (2.1)--(2.2) is a correct
   measure-valued branching transport in EGF space.
2. **Signed coefficient extraction.** The boundary is (3.5), a bivariate
   Cauchy coefficient of an indefinite double-resolvent kernel.
3. **Residue cross-variation.** At a pole $\mu$ of $q$, the cross-variation is
   (3.3), i.e. $-\mu$ times the renewal residue mass $G_q(\mu)$.
4. **Missing ingredient.** A proof of (G') requires a domination theorem for the
   signed residue masses of the chapter renewal at negative poles of $O_{n-d}$.
   Equivalently, a positive martingale after a **signed** change of measure, or
   a de Branges/Krein-type structural identity, not a positive chapter coupling.

No such domination theorem is proved here, and the constructions above show why
the natural positive-coupling routes cannot supply it.

---

## 8. Files and verification

- `NOTES.md` — this report.
- `verify_coupling.py` — exact symbolic checks: CB/Möbius flow identity
  (2.2)/(2.4), Wronskian B1 identity, residue identity (3.4)--(3.5),
  covariance identity (1.6) for $t>0$, sampled boundary
  $W(H_n,O_{n-d})$ real-root absence, and $\gamma>0$ on $(0,z_B)$ for
  $j=2,\ldots,9$. Run `py -3 verify_coupling.py` from this directory.

No external packages are required beyond `sympy` and Python stdlib. No builds,
no large caches.

## 9. Final conclusion

The measure-valued branching coupling suggested by
`proof-chapter-renewal/NOTES.md` §5 is **constructed**:

1. **CB flow coupling.** The quadratic inhomogeneous CB Riccati flow (2.1)--(2.3)
   is a positive continuous-state branching coupling for the Möbius flow
   $\phi_z$, with affine chapter observable (2.5).
2. **Spine coupling.** The renewal spine (3.1) is a positive discrete
   martingale decomposition of $H_n$ for $t>0$.
3. **Pole cross-variation.** The signed residue mass (3.2) gives the exact pole
   cross-variation $\mathcal X_{n,d}(\mu)=-\mu G_q(\mu)$ in (3.5)--(3.6),
   which is equivalent to the boundary Pick inequality (G).

The **proof is not completed** because positivity of $\mathcal X_{n,d}(\mu)$
remains open. The obstruction is precise:

- positive CB/flow positivity does not survive signed coefficient extraction
  (Theorem 4.1, Corollary 4.2);
- the reciprocal kernel $1/g=B/C$ is not globally monotone/Stieltjes
  (Theorem 5.3);
- the residue mass $\Xi_{n,d,\mu}$ is not termwise positive (Section 6, G').

Therefore the requested construction exists, but it yields a rigorous signed
cross-variation identity and a precise obstruction, not a proof of
$W(H_n,O_{n-d})>0$.