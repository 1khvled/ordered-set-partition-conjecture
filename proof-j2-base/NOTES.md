# `j=2` base case: exact reduction, limiting sequence, and obstructions

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-j2-base` was created/edited.

## Outcome

**No unconditional all-`n` proof of the `j=2` base case was obtained.**

The new rigorous content is:

1. An exact two-variable deformation of the consecutive residue, reducing the
   target at `x=1` to a root-velocity inequality.
2. An all-index calculation of the `x -> infinity` limit.  After normalization
   it is the classical Fibonacci-polynomial recurrence, whose roots are simple,
   positive, and consecutively interlacing.
3. A precise continuation criterion: two discriminant/resultant positivity
   statements for the deformed family would finish the proof.
4. An associated-Stirling ratio lemma, proved from the strong log-concavity of
   the no-singleton Bell polynomials.
5. Exact obstructions to several natural associated-Stirling/Fubini,
   falling-basis, and quadratic-block shortcuts.

The original targets remain

\[
\frac{O_{n+1}(1,\mu)}{O_n'(1,\mu)}<0
\quad\text{at every root }\mu\text{ of }O_n(1,t),
\]

and

\[
W(O_n,O_{n-2})(t)>0\qquad(t\le 0).
\]

No finite check below is used as an all-index proof.

---

## 1. The `j=2` family and a correct exact residue identity

Put

\[
C(z)=e^z-1-z,\qquad
F(z,t)=\frac1{1-z-tC(z)},\qquad
O_n(t)=n![z^n]F(z,t).
\]

Let

\[
U_n(t)=n![z^n]F(z,t)^2
      =\sum_{r=0}^n\binom nr O_r(t)O_{n-r}(t).
\]

Two elementary identities are

\[
F_t=CF^2,
\qquad
F_z=\bigl(2+(t-1)z\bigr)F^2-F.
\]

Since `tC=(1-z)-F^{-1}`, the first identity gives

\[
tO_n'=U_n-nU_{n-1}-O_n. \tag{1}
\]

The second gives

\[
O_{n+1}=2U_n+(t-1)nU_{n-1}-O_n. \tag{2}
\]

Eliminating `U_n` between (1) and (2) yields the useful exact formula

\[
\boxed{
O_{n+1}=O_n+2tO_n'+n(1+t)U_{n-1}.
} \tag{3}
\]

Consequently, at a root `mu` of `O_n`,

\[
\boxed{
\frac{O_{n+1}(\mu)}{O_n'(\mu)}
=2\mu+n(1+\mu)\frac{U_{n-1}(\mu)}{O_n'(\mu)}.
} \tag{4}
\]

This is a `j=2`-only reduction.  It is exact and does not use the generic
cross-row or coefficient-minor arguments recorded elsewhere in the repository.

A sufficient condition is

\[
\frac{U_{n-1}(\mu)}{O_n'(\mu)}>0,
\]

but (4) shows that a quantitative upper bound is needed when `-1<mu<0`.
Positivity alone does not close the target.

---

## 2. Deforming the singleton weight

Introduce a second parameter `x`:

\[
F(x,z,t)=\frac1{1-xz-tC(z)},\qquad
O_n(x,t)=n![z^n]F(x,z,t).
\]

The original family is `x=1`.  Because

\[
\partial_xF=zF^2,
\]

we have

\[
\partial_xO_n(x,t)=nU_{n-1}(x,t). \tag{5}
\]

Let `mu_n(x)` be a simple root branch of `O_n(x,t)`.  Implicit differentiation
and (5) give

\[
\mu_n'(x)=-\frac{nU_{n-1}(x,\mu_n(x))}
                   {O_n'(x,\mu_n(x))}.
\]

Substituting this into (4) gives the exact root-velocity form

\[
\boxed{
\frac{O_{n+1}(1,\mu)}{O_n'(1,\mu)}
=2\mu-(1+\mu)\mu_n'(1).
} \tag{6}
\]

Thus the consecutive target is equivalent to

\[
\boxed{
(1+\mu)\mu_n'(1)>2\mu
\quad\text{for every root }\mu\text{ of }O_n(1,t).
} \tag{7}
\]

This is a concrete `j=2`-specific compatible-sequence question: all roots must
move left under the singleton-weight deformation at precisely the rate in
(7).

---

## 3. The all-index limit `x -> infinity`

Scale

\[
t=x^2u.
\]

From the explicit coefficient formula one obtains

\[
\lim_{x\to\infty}x^{-n}O_n(x,x^2u)=\Phi_n(u),
\]

where

\[
\Phi_n(u)=
\sum_{k=0}^{\lfloor n/2\rfloor}
\frac{n!(n-k)!}{(n-2k)!\,2^k k!}\,u^k. \tag{8}
\]

Writing `phi_n=Phi_n/n!`, the first-block recurrence gives

\[
\phi_{n+1}(u)=\phi_n(u)+\frac u2\phi_{n-1}(u),
\qquad \phi_0=\phi_1=1. \tag{9}
\]

With `v=-u/2`, define

\[
P_n(v)=\phi_n(-2v).
\]

Then

\[
P_{n+1}(v)=P_n(v)-vP_{n-1}(v),
\qquad P_0=P_1=1,
\]

or explicitly

\[
P_n(v)=\sum_{j=0}^{\lfloor n/2\rfloor}
(-1)^j\binom{n-j}{j}v^j. \tag{10}
\]

These are the classical Fibonacci polynomials.  Their roots are simple and
positive, and their consecutive root pattern is:

* `P_{2m}` has `m` roots and surrounds the `m-1` roots of `P_{2m-1}`;
* `P_{2m+1}` has `m` roots, with its `i`th root to the left of the `i`th root
  of `P_{2m}`.

This follows by induction from `P_{n+1}=P_n-vP_{n-1}`: at a root of `P_n`,
`P_{n+1}=-vP_{n-1}`, and the induction hypothesis supplies the alternating
signs; Rolle's theorem supplies the roots between consecutive roots, while the
signs at `0` and `+infinity` supply the remaining root(s).

Since `u=-2v`, every `Phi_n` has simple negative roots.  The interlacing
orientation in `u` is exactly the one required by the consecutive sign
criterion, so

\[
\frac{\Phi_{n+1}(u)}{\Phi_n'(u)}<0
\quad\text{at every root }u\text{ of }\Phi_n. \tag{11}
\]

Moreover,

\[
\frac{O_{n+1}(x,x^2u)}{O_n'(x,x^2u)}
\sim x^3\frac{\Phi_{n+1}(u)}{\Phi_n'(u)},
\]

so the target residue is negative for all sufficiently large `x`.

---

## 4. A precise continuation criterion

The limit result suggests a finite-`x` continuation argument.  It is enough to
prove the following two statements for the deformed family.

### (D) No root collision for every positive `x`

For every `n` and `x>0`, the discriminant of `O_n(x,t)` in `t` is nonzero
(indeed, it is observed to be positive).  Together with the simple real limit
at `x -> infinity`, this keeps every root real and simple for all `x>0`.
Once the roots are real, they are negative because the coefficients are
positive and `O_n(x,0)>0`.

### (R) No consecutive common root for every positive `x`

For every `n` and `x>0`,

\[
\operatorname{Res}_t\bigl(O_n(x,t),O_{n+1}(x,t)\bigr)\ne0.
\]

If (D) and (R) hold, the roots `mu_n(x)` continue as simple real branches for
all `x>0`.  By (6), the residue can change sign only if it becomes zero, i.e.
only if `O_n` and `O_{n+1}` acquire a common root.  Condition (R) excludes
this.  Its sign is negative at `x -> infinity` by (11), hence negative at
`x=1`.

Thus (D)+(R) are a complete all-index route to the consecutive target.

### Exact finite evidence for (D) and (R)

For example,

\[
\operatorname{Res}_t(O_2,O_3)=-2x^2(3x+1),
\]

\[
\operatorname{Res}_t(O_3,O_4)
=-6x^3(36x^3+36x^2+10x+1),
\]

and

\[
\operatorname{Res}_t(O_4,O_5)
=-48x^4(10800x^6+21600x^5+16200x^4+5920x^3
+1155x^2+130x+7).
\]

Also,

\[
\operatorname{Disc}_t O_4
=(2x+1)(6x+1)(60x^2+8x+1),
\]

and

\[
\operatorname{Disc}_t O_5
=14400x^6+19200x^5+8400x^4+1680x^3+220x^2+20x+1.
\]

The exact script `verify.py` checks these and the analogous computations
through `n=6`.  All observed residual factors have positive coefficients, but
this finite observation is not a proof of (D) or (R).

---

## 5. Associated-Stirling structure

Let

\[
A_{m,k}={m\brace k}_{\ge2}
\]

be the number of partitions of `[m]` into `k` blocks, all of size at least two.
Then

\[
A_{m,k}=kA_{m-1,k}+(m-1)A_{m-2,k-1}. \tag{12}
\]

The monomial coefficients of the original row are

\[
[t^k]O_n(1,t)
=n!\sum_{m=2k}^n
\frac{(n-m+1)_k}{m!}A_{m,k}. \tag{13}
\]

### 5.1 A proved associated-Stirling ratio lemma

Bóna--Mező's theorem says that the row polynomial
`sum_k A_{m,k}t^k` has simple nonpositive roots.  Its Newton inequalities
give, in particular,

\[
A_{m,k}^2\ge \frac{k+1}{k}A_{m,k-1}A_{m,k+1}. \tag{14}
\]

Define

\[
\Delta_{m,k}=A_{m-1,k}A_{m,k+1}
             -A_{m-1,k+1}A_{m,k}.
\]

Using (12), direct expansion gives

\[
\Delta_{m,k}
=A_{m-1,k}A_{m-1,k+1}+(m-1)E_{m,k}, \tag{15}
\]

where

\[
\begin{aligned}
E_{m,k}
={}&A_{m-2,k}A_{m-1,k}
   -A_{m-2,k-1}A_{m-1,k+1}\\
={}&kA_{m-2,k}^2-(k+1)A_{m-2,k-1}A_{m-2,k+1}
   +(m-2)\Delta_{m-2,k-1}.
\end{aligned} \tag{16}
\]

The first two terms in (16) are nonnegative by (14).  Induction on `m` in
(15)--(16) therefore proves

\[
\boxed{
\frac{A_{m,k+1}}{A_{m,k}}
\quad\text{is increasing in }m
} \tag{17}
\]

on the active range.  This is an all-index associated-Stirling result.

### 5.2 Falling-factorial coefficients and an exact product formula

Write

\[
O_n(1,t)=\sum_i c_{n,i}(t)_i,
\qquad (t)_i=t(t-1)\cdots(t-i+1).
\]

Then

\[
c_{n,i}=\sum_{k\ge i}{k\brace i}[t^k]O_n(1,t). \tag{18}
\]

A finite-difference identity gives the more useful exact formula

\[
\boxed{
c_{n,i}
=n![z^n]\frac{C(z)^i}
 {\prod_{r=0}^i(1-z-rC(z))}.
} \tag{19}
\]

Indeed, for `f(j)=1/(1-z-jC)`,

\[
\frac1{i!}\Delta_j^i f(0)
=\frac{C^i}{\prod_{r=0}^i(1-z-rC)}.
\]

The natural compatible-sequence guess is that
`c_{n+1,i}/c_{n,i}` increases with `i`.  Exact checks support it, but no
all-index proof was obtained.

### 5.3 A precise missing associated-Stirling inequality

A stronger inequality that would imply useful TP2 properties for the matrix
in (13) is

\[
\boxed{
 kA_{m+1,k+1}A_{m,k}
 \ge (k+1)A_{m+1,k}A_{m,k+1}.
} \tag{H}
\]

It is equivalent to

\[
\frac{A_{m+1,k+1}/A_{m+1,k}}
     {A_{m,k+1}/A_{m,k}}
\ge \frac{k+1}{k}.
\]

Exact checks pass through `m<=300`, `k<=100`, but (H) is not proved here.
The proved lemma (17) is its unweighted version and does not imply (H).

Nor does ordinary Bóna--Mező interlacing imply (H).  An exact counterexample
to that inference is

\[
P=(t+1)(t+2)(t+3),
\qquad
Q=(t+\tfrac32)(t+\tfrac52)(t+4).
\]

The roots of `Q` alternate left of those of `P`, but with
`p_i=[t^i]P`, `q_i=[t^i]Q`,

\[
q_2p_1-2q_1p_2=-149<0.
\]

Thus invoking the associated-Bell interlacing theorem alone cannot supply
(H).

### 5.4 Why simple falling-basis ratio monotonicity is insufficient

Even if `c_{n+1,i}/c_{n,i}` were proved increasing, that condition by itself
does not force the Bezoutian to be positive definite.  For example,

\[
P=10+11t+t^2,\qquad Q=10+22t+100t^2
\]

have falling-factorial coefficient vectors

\[
(10,12,1),\qquad (10,122,100),
\]

whose componentwise ratios are strictly increasing, while

\[
\operatorname{Disc}_t(Q)=22^2-4\cdot10\cdot100=-3516<0.
\]

So a proof based only on the natural falling-basis ratio condition is
blocked unless an additional banded/minor condition is supplied.

---

## 6. The natural no-singleton Fubini intermediate does not directly interlace

The standard intermediate

\[
G(z,t)=\frac1{1-tC(z)}
\]

satisfies

\[
F(x,z,t)=\frac{G(z,t)}{1-xzG(z,t)}.
\]

The known ordered-lift theorem applies to `e^zG`, not to `G` or to the
required geometric transform.

Moreover, `G` is not directly compatible with the target row.  At `n=4`,

\[
G_4=t+6t^2
\]

has roots `-1/6,0`, whereas

\[
O_4=24+45t+6t^2
\]

has roots

\[
\frac{-45\pm\sqrt{1449}}{12},
\]

both below `-1/2`.  The two root sets do not interlace in either same-degree
orientation.  This blocks a direct preservation-of-interlacing argument from
`G_n` to `O_n`.

---

## 7. Exact obstruction to a simple quadratic-block Wronskian certificate

For lag two, put

\[
W_n(t)=W(O_n,O_{n-2})(t)=\sum_k w_{n,k}t^k.
\]

A tempting negative-axis certificate is to split `W_n(-s)` into disjoint
positive quadratics

\[
w_k-w_{k+1}s+w_{k+2}s^2.
\]

This already fails by an exact coefficient inequality at `n=8`.  Here

\[
(w_{8,3},w_{8,4},w_{8,5})
=(323112180,106099140,7509600),
\]

but

\[
106099140^2
=11257027508739600
>
9705772907712000
=4\cdot323112180\cdot7509600.
\]

Thus that quadratic has positive real roots.  This does not contradict WRS;
it only blocks this simple fixed grouping certificate.

---

## 8. What remains to be proved

The cleanest remaining all-index package is (D)+(R) from Section 4.  It would
imply the full consecutive target by continuation from the Fibonacci limit.

Alternatively, a proof must establish the exact root-velocity inequality (7),
or a genuinely stronger compatible-sequence statement than the insufficient
falling-basis ratio condition of Section 5.4.  The associated-Stirling
inequality (H) is one concrete possible bottleneck, but it is not known and is
not a consequence of the standard associated-Bell interlacing theorem.

The lag-two target

\[
W(O_n,O_{n-2})(t)>0\qquad(t\le0)
\]

was not proved.  The exact identity (3) concerns consecutive rows; it does not
by itself control the lag-two residues.

---

## 9. Files and reproduction

* `NOTES.md` -- this report.
* `verify.py` -- exact integer/rational checks for (1)--(5), the associated
  Stirling lemma on a finite range, the small resultants/discriminants, and the
  exact obstructions above.
* `explore.py` -- small numerical probes used during the search; not a proof.

Run from this directory:

```powershell
py -3.14 verify.py
```

No packages were installed, no builds were run, and all created files are
small.
