# Associated-Stirling ratio inequality (H) — algebraic reduction and obstruction

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-associated-stirling` was created/edited.

## Target

Let
$$
A_{m,k} := {m \brace k}_{\ge 2}
$$
be the 2-associated Stirling number of the second kind: the number of partitions
of $[m]$ into $k$ blocks, every block of size at least $2$.  The active range is
$0\le k\le \lfloor m/2\rfloor$, with $A_{0,0}=1$ and $A_{m,0}=0$ for $m>0$.
The recurrence is
$$
A_{m,k}=kA_{m-1,k}+(m-1)A_{m-2,k-1},\qquad m\ge1,k\ge1, \tag{1}
$$
with the boundary convention $A_{m,k}=0$ when $2k>m$ or $k<0$.

The inequality identified in `proof-j2-base/NOTES.md` §5.3 is
$$
\boxed{
D_{m,k}:=k\,A_{m+1,k+1}A_{m,k}-(k+1)A_{m+1,k}A_{m,k+1}\ge 0.
} \tag{H}
$$
Equivalently, with $\rho(m,k):=A_{m,k+1}/(kA_{m,k})$ on the active range,
$$
\rho(m+1,k)\ge \rho(m,k). \tag{H'}
$$

## Outcome

**No unconditional all-index proof of (H) was obtained.**

Rigorous content obtained:

1. An exact algebraic reduction of (H) to a single-row strong log-concavity
   term plus a nonnegative forcing term.
2. A two-step recurrence for the reduced quantity $E'_{m,k}$ whose forcing term
   is numerically positive on a large range but whose SLC component is not
   implied by Bóna–Mező real-rootedness alone.
3. Exact finite verification of (H), the reduction identities, and the forcing
   term through $m\le 300$ (and larger spot checks reported in script output).

No counterexample to (H) was found.

---

## 1. Proven identities

### 1.1 The one-step reduction

Substituting (1) in the form
$$
A_{m+1,k+1}=(k+1)A_{m,k+1}+mA_{m-1,k},\qquad
A_{m+1,k}=kA_{m,k}+mA_{m-1,k-1}
$$
into (H) gives
$$
\begin{aligned}
D_{m,k}
&=k\bigl((k+1)A_{m,k+1}+mA_{m-1,k}\bigr)A_{m,k}
 -(k+1)\bigl(kA_{m,k}+mA_{m-1,k-1}\bigr)A_{m,k+1}\\
&=m\bigl(kA_{m-1,k}A_{m,k}-(k+1)A_{m-1,k-1}A_{m,k+1}\bigr).
\end{aligned}
$$
Define
$$
E'_{m,k}:=kA_{m-1,k}A_{m,k}-(k+1)A_{m-1,k-1}A_{m,k+1}. \tag{2}
$$
Thus
$$
\boxed{D_{m,k}=m\,E'_{m,k},\qquad m\ge1,} \tag{3}
$$
and (H) is equivalent to $E'_{m,k}\ge0$ on its active range.

This identity is algebraic and exact (verified symbolically in `verify.py`).

### 1.2 The two-step recurrence for $E'$

Expand $A_{m,k}$ and $A_{m,k+1}$ in (2) using (1):
$$
A_{m,k}=kA_{m-1,k}+(m-1)A_{m-2,k-1},\quad
A_{m,k+1}=(k+1)A_{m-1,k+1}+(m-1)A_{m-2,k}.
$$
A direct expansion yields
$$
\begin{aligned}
E'_{m,k}
={}&k^2A_{m-1,k}^2-(k+1)^2A_{m-1,k-1}A_{m-1,k+1}\\
&+(m-1)\bigl(kA_{m-1,k}A_{m-2,k-1}-(k+1)A_{m-1,k-1}A_{m-2,k}\bigr).
\end{aligned} \tag{4}
$$
Introduce the single-row strong log-concavity term
$$
SLC(r,k):=k^2A_{r,k}^2-(k+1)^2A_{r,k-1}A_{r,k+1}, \tag{5}
$$
the unweighted ratio difference
$$
\Delta(r,k):=A_{r-1,k}A_{r,k+1}-A_{r-1,k+1}A_{r,k}, \tag{6}
$$
and the weighted cross difference appearing in the bracket of (4),
$$
G_{m,k}:=kA_{m-1,k}A_{m-2,k-1}-(k+1)A_{m-1,k-1}A_{m-2,k}.
$$
The identity
$$
G_{m,k}=D_{m-2,k-1}+\Delta(m-1,k-1) \tag{7}
$$
is exact: with $D_{p,j}=jA_{p+1,j+1}A_{p,j}-(j+1)A_{p+1,j}A_{p,j+1}$ one gets
$D_{m-2,k-1}=(k-1)A_{m-1,k}A_{m-2,k-1}-kA_{m-1,k-1}A_{m-2,k}$, and adding
$\Delta(m-1,k-1)=A_{m-2,k-1}A_{m-1,k}-A_{m-2,k}A_{m-1,k-1}$ gives $G_{m,k}$.

Using (3) in the form $D_{m-2,k-1}=(m-2)E'_{m-2,k-1}$, we obtain the
two-step recurrence
$$
\boxed{
E'_{m,k}=SLC(m-1,k)+(m-1)\Delta(m-1,k-1)+(m-1)(m-2)E'_{m-2,k-1}.
} \tag{8}
$$

The forcing term is
$$
F(m,k):=SLC(m-1,k)+(m-1)\Delta(m-1,k-1). \tag{9}
$$
Thus $E'_{m,k}=F(m,k)+(m-1)(m-2)E'_{m-2,k-1}$.

All identities (3)–(9) are verified symbolically/exactly in `verify.py`.

### 1.3 Known nonnegative ingredients

* **$\Delta(r,k)\ge0$** is proved in `proof-j2-base/NOTES.md` §5.1, as a
  consequence of Bóna–Mező real-rootedness (Newton inequality
  $A_{r,k}^2\ge \frac{k+1}{k}A_{r,k-1}A_{r,k+1}$) and induction via the
  same recurrence (1).
* **$E'_{m,k}\ge0$ for base cases** is checked exactly by the script
  (and follows from the small-$m$ values of the associated Stirling table).

Therefore (H) follows for all $m$ by induction on $m$ in steps of 2 once
$F(m,k)\ge0$ is known on the active range.

---

## 2. The obstruction

### 2.1 SLC is stronger than Newton and not a real-rootedness consequence

Bóna–Mező prove that $B_r(x)=\sum_k A_{r,k}x^k$ is real-rooted of degree
$d=\lfloor r/2\rfloor$. Newton's inequality gives
$$
A_{r,k}^2\ge \frac{k+1}{k}\frac{d-k+1}{d-k}A_{r,k-1}A_{r,k+1}. \tag{10}
$$
This implies SLC only when $\frac{d-k+1}{d-k}\ge\frac{k+1}{k}$, i.e.
$2k\ge d$ (the upper half). In the lower half $2k<d$, Newton gives a weaker
factor and cannot imply SLC.

Moreover SLC is false for general real-rooted polynomials: for
$(1+x)^d$, $k=1$, one has $c_1^2/(c_0c_2)=2d/(d-1)<4$ for $d\ge3$, while
SLC would require $\ge4$. Hence SLC is a special property of the associated
Stirling sequence, not a formal consequence of real-rootedness.

### 2.2 Direct SLC induction does not close

Expanding $SLC(m,k)$ by (1) with $M=m-1$,
$A=A_{m-1,k},B=A_{m-1,k-1},C=A_{m-1,k+1},P=A_{m-2,k-1},Q=A_{m-2,k-2},R=A_{m-2,k}$,
gives
$$
SLC(m,k)=M^2(k^2P^2-(k+1)^2QR)+M(\cdots)+
 (k^4A^2-(k^2-1)(k+1)^2BC).
$$
Neither the $M^2$ coefficient nor the $M^0$ coefficient is the same SLC term
required by a direct induction; the $M^0$ coefficient is even more demanding
than SLC. Therefore the straightforward recurrence induction for SLC does not
close.

### 2.3 Status

* (H) is **not proved** here.
* SLC is **not proved** here.
* $F(m,k)\ge0$ is **not proved** here.
* Numerical evidence (no counterexamples) was found for (H), SLC, and $F$ on
  large finite ranges; see `verify.py`.
* The exact reduction (8) means: a proof of $SLC(m-1,k)\ge0$, or even the
  weaker $F(m,k)\ge0$, on the active range would immediately prove (H).

---

## 3. Verification

`verify.py` performs exact integer/rational checks for:
- recurrence (1) table construction;
- identity (3) $D=mE'$;
- identity (8) two-step recurrence for $E'$;
- (H) on a finite range;
- SLC on a finite range;
- $F$ on a finite range;
- the binomial real-rootedness counterexample to "Newton => SLC".

Run from this directory:

```powershell
py -3.14 verify.py
```
