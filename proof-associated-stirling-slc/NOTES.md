# Associated-Stirling SLC: refined recurrence, G_r real-rootedness, and obstruction

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-associated-stirling-slc` was created/edited.

## Target

For 2-associated Stirling numbers $A_{m,k}=\{m\brace k\}_{\ge2}$, prove
$$
SLC(r,k):=k^2A_{r,k}^2-(k+1)^2A_{r,k-1}A_{r,k+1}\ge0
$$
on the active range $1\le k\le \lfloor r/2\rfloor-1$, or prove the weaker
$F(m,k):=SLC(m-1,k)+(m-1)\Delta(m-1,k-1)\ge0$ which implies (H) by identity (8)
of `proof-associated-stirling/NOTES.md`.

Equivalently, $SLC(r,k)\ge0$ is log-concavity of $c_{r,k}:=(k!)^2A_{r,k}$.

## Outcome

**No unconditional all-index proof of SLC or $F\ge0$ was obtained.**

New rigorous content:

1. A sufficient real-rootedness condition: if $G_r(x):=\sum_{k=0}^{\lfloor r/2\rfloor}(k!)^2A_{r,k}x^k$ is real-rooted, Newton gives SLC. $G_r$ is real-rooted for $r\le40$ (exact/MP check, finite evidence).
2. A corrected/explicit SLC recurrence expansion showing the $M^2$ and $M^0$ coefficients are *implied* by SLC on previous rows; the remaining obstruction is a single cross term $C_1$.
3. An exact obstruction to the naive row-log-concavity induction: the recurrence $c_{m,k}=k c_{m-1,k}+(m-1)k^2c_{m-2,k-1}$ does not preserve log-concavity by itself.

No counterexample to SLC, $F$, or $G_r$ real-rootedness was found.

---

## 1. The G_r sufficient condition

Let $c_{r,k}=(k!)^2A_{r,k}$ and $G_r(x)=\sum_{k=0}^{d_r}c_{r,k}x^k$, $d_r=\lfloor r/2\rfloor$.
SLC is $c_{r,k}^2\ge c_{r,k-1}c_{r,k+1}$.

If $G_r$ is real-rooted, Newton's inequality gives
$$
 c_{r,k}^2\ge \frac{k+1}{k}\frac{d_r-k+1}{d_r-k}c_{r,k-1}c_{r,k+1}\ge c_{r,k-1}c_{r,k+1},
$$
hence SLC. This is strictly stronger than SLC.

Using $A_{r,k}=\frac{r!}{k!}[z^r](e^z-1-z)^k$, one has
$$
 G_r(x)=r![z^r]\sum_{k\ge0}k!(x(e^z-1-z))^k
$$
as a finite formal coefficient extraction (finite for each $r$ since $e^z-1-z=O(z^2)$).
The recurrence $A_{m,k}=kA_{m-1,k}+(m-1)A_{m-2,k-1}$ gives
$$
 c_{m,k}=k c_{m-1,k}+(m-1)k^2 c_{m-2,k-1},
$$
and hence the differential recurrence
$$
 G_m(x)=xG_{m-1}'(x)+(m-1)x(1+xD)^2G_{m-2}(x),\qquad G_0=1,G_1=0,\ D=x\partial_x.
$$

Finite evidence: $G_r$ is real-rooted for $2\le r\le40$ (mpmath, 80 dps; `verify.py`).
Also $G_r$ has nonpositive roots and a factor $x$ for $r\ge2$.

**Gap:** proving real-rootedness of $G_r$ for all $r$ is open here. The two-step
operator $(1+xD)^2$ is not a standard real-rootedness-preserving composition in
this form, and a direct interlacing induction was not closed.

---

## 2. Refined SLC recurrence expansion

Let $M=m-1$, $A=A_{m-1,k},B=A_{m-1,k-1},C=A_{m-1,k+1}$, $P=A_{m-2,k-1},Q=A_{m-2,k-2},R=A_{m-2,k}$.
Using $A_{m,k}=kA+MP$, $A_{m,k-1}=(k-1)B+MQ$, $A_{m,k+1}=(k+1)C+MR$ gives
$$
 SLC(m,k)=M^2C_2(m,k)+M C_1(m,k)+C_0(m,k),
$$
where
$$
\begin{aligned}
C_2&=k^2P^2-(k+1)^2QR,\\
C_1&=2k^3AP-(k+1)^2\bigl((k-1)BR+(k+1)CQ\bigr),\\
C_0&=k^4A^2-(k^2-1)(k+1)^2BC.
\end{aligned}
$$

### 2.1 C0 and C2 are controlled by SLC

Exact identities:
$$
C_0=k^2\,SLC(m-1,k)+(k+1)^2BC,
$$
$$
C_2=SLC(m-2,k-1)+(2k-1)P^2-(2k+1)QR.
$$
Moreover $SLC(m-2,k-1)=(k-1)^2P^2-k^2QR\ge0$ implies $P^2/QR\ge k^2/(k-1)^2>(k+1)^2/k^2$, hence $C_2\ge0$.
Also $SLC(m-1,k)\ge0$ implies $C_0\ge0$.

Thus, under induction hypotheses $SLC(m-1,k)\ge0$ and $SLC(m-2,k-1)\ge0$, the
quadratic in $M$ has nonnegative end coefficients.

**Correction to note 2.2 in `proof-associated-stirling/NOTES.md`:** that note
stated $C_0$ is "even more demanding than SLC"; the identity above shows the
opposite: $C_0\ge k^2SLC(m-1,k)$.

### 2.2 The remaining cross term

The unresolved term is $C_1$. Finite checks find no $C_1<0$ for $m\le120$
(worst ratio $C_1/[(k+1)^2((k-1)BR+(k+1)CQ)]\approx0.150$ at $(m,k)=(120,36)$;
`verify.py`).

However $C_1\ge0$ is not implied by SLC/Delta on previous rows alone in an
obvious ratio form: writing $u_k=A_{m-1,k}/A_{m-1,k-1}$, $v_k=A_{m-2,k}/A_{m-2,k-1}$,
$$
\frac{C_1}{AP}=2k^3-(k+1)^2\left((k-1)\frac{v_k}{u_k}+(k+1)\frac{u_{k+1}}{v_{k-1}}\right),
$$
and the available bounds ($v_k/u_k\le1$, $u_{k+1}/u_k\le k^2/(k+1)^2$, $v_k/v_{k-1}\le (k-1)^2/k^2$)
do not control the $u_k/v_k$ factor inside $u_{k+1}/v_{k-1}$.

Therefore the SLC recurrence induction is reduced to the single inequality
$C_1(m,k)\ge0$ on the active range, but $C_1$ is not closed by the present
auxiliary invariants.

---

## 3. Obstruction: naive LC propagation fails

The row recurrence $c_{m,k}=k c_{m-1,k}+(m-1)k^2c_{m-2,k-1}$ does not preserve
log-concavity of arbitrary nonnegative sequences. An explicit counterexample is
in `verify.py` (`naive LC propagation: counterexample confirmed`). Hence any
induction proving SLC must use additional structure of $A_{r,k}$ beyond row LC.

This rules out the direct "LC in, LC out" closure for the natural $c_{r,k}$
recurrence.

---

## 4. Finite verification (not a proof)

`verify.py` checks, for $m\le120$ unless noted:

- SLC: all checked entries nonnegative (worst ratio $\approx0.240$ at $(80,24)$).
- SLC expansion identity $SLC(m,k)=M^2C_2+MC_1+C_0$.
- $C_0=k^2SLC(m-1,k)+(k+1)^2BC$ and $C_2$ implied by prior SLC.
- $C_1>0$ (worst ratio $\approx0.150$ at $(120,36)$).
- $G_r$ real-rooted for $r\le40$ (mpmath 80 dps).
- Naive LC propagation counterexample.

Run:

```powershell
py -3.14 verify.py
```

---

## 5. Precise remaining gap

A proof of any of the following would close SLC (and hence, via (8) in
`proof-associated-stirling/NOTES.md`, close (H)):

1. **$G_r$ real-rooted for all $r$** (sufficient, stronger).
2. **$C_1(m,k)\ge0$ for all active $m,k$** under the SLC induction (sufficient by 2.1).
3. **$F(m,k)\ge0$** directly (weaker than SLC, sufficient for (H)).
4. A **refined auxiliary invariant** stronger than SLC but preserved by the
   $(k, (m-1)k^2)$ recurrence, or a PF∞/combinatorial argument for
   $c_{r,k}=(k!)^2A_{r,k}$.

The all-index statement remains open. The strongest new reduction here is:
$SLC(m-1,\cdot),SLC(m-2,\cdot)\ge0 \implies SLC(m,k)\ge0$ provided $C_1(m,k)\ge0$,
with $C_0,C_2$ rigorously controlled and $C_1$ an explicit single cross-term obstruction.
