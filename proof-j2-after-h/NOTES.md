# j=2 base case after (H): does (H) close (7) or the TP2/Bezoutian route?

Date: 2026-11-28
Scope: only `ordered-big-blocks/proof-j2-after-h` created/edited.
No installs/builds. No other agents/models.

## Outcome

**No closure obtained.** Inequality (H) for 2-associated Stirling numbers,
now proved unconditionally in `proof-u-operator/NOTES.md` §7, together with
the proved unweighted ratio lemma / associated-Stirling ratio monotonicity
(17) from `proof-j2-base/NOTES.md` §5.1, does **not** by itself close:

1. the j=2 root-velocity inequality (7) from `proof-j2-base/NOTES.md` §2, nor
2. the TP2/Bezoutian positivity needed for the j=2 falling-factorial route.

The exact remaining gap is identified below, with a feasible falsification
check (finite evidence only, not a proof).

---

## 1. Inputs used (all previously proved / stated)

### (H) — proved, `proof-u-operator/NOTES.md` §7

For $A_{m,k}={m\brace k}_{\ge2}$, $m=r-1$, $k=j+1$:
$$
 k A_{m+1,k+1}A_{m,k}\ge (k+1)A_{m+1,k}A_{m,k+1}. \tag{H}
$$
Equivalently
$$
 \frac{A_{m+1,k+1}/A_{m+1,k}}{A_{m,k+1}/A_{m,k}}\ge \frac{k+1}{k}. \tag{H'}
$$
Proved by a coupled induction on the row polynomials $G_r(x)=E_r(x)/x=\sum_{k\ge1}k!A_{r,k}x^{k-1}$:
$G_{r-1}\preceq G_r$ and $\mathcal U G_{r-1}\preceq G_r$, then the coefficient-ratio
inequality for interlacing nonpositive-root polynomials gives (H).
The proof also states: "This proves (H)... it does **not** by itself prove the
full $j=2$ root-velocity or discriminant/resultant statements."

### (17) — proved, `proof-j2-base/NOTES.md` §5.1

$$
 \frac{A_{m,k+1}}{A_{m,k}}\quad\text{is increasing in }m
 \tag{17}
$$
on the active range. This is the unweighted version of (H): it follows from
Bóna–Mező Newton inequality (14) and the recurrence (15)–(16). Note (H) is
strictly stronger than (17) since $(k+1)/k>1$.

### (7) — target, `proof-j2-base/NOTES.md` §2

With $\mu_n(x)$ a simple root branch of $O_n(x,t)$,
$$
 \frac{O_{n+1}(1,\mu)}{O_n'(1,\mu)}=2\mu-(1+\mu)\mu_n'(1), \tag{6}
$$
and the consecutive target is equivalent to
$$
 (1+\mu)\mu_n'(1)>2\mu
 \quad\text{for every root }\mu\text{ of }O_n(1,t). \tag{7}
$$

### Falling-factorial / TP2 / Bezoutian route (j2-base §5.2–5.4)

Monomial coefficients (13):
$$
 d_{n,k}:=[t^k]O_n(1,t)=n!\sum_{m=2k}^{n}\frac{(n-m+1)_k}{m!}A_{m,k}. \tag{13}
$$
Falling-factorial coefficients (18)–(19):
$$
 O_n(1,t)=\sum_i c_{n,i}(t)_i,\qquad
 c_{n,i}=\sum_{k\ge i}{k\brace i}d_{n,k}. \tag{18}
$$
Natural guess: $c_{n+1,i}/c_{n,i}$ increasing in $i$. 5.4 obstruction:
ratio monotonicity alone does not force Bezoutian positive definiteness
(counterexample $P=10+11t+t^2$, $Q=10+22t+100t^2$, increasing falling-basis
ratios but $\operatorname{Disc}_t(Q)=-3516<0$). Thus an additional
banded/minor condition beyond ratio monotonicity is required.

---

## 2. What (H)+(17) give, and what they do not

### 2.1 Immediate consequences of (H)+(17)

- (17) is unweighted TP2 of $A_{m,k}$ in adjacent $(m,k)$:
  $A_{m+1,k+1}A_{m,k}\ge A_{m+1,k}A_{m,k+1}$.
- (H) is stronger: multiplies the right side by $(k+1)/k>1$.
- Both are coefficient-level inequalities on the associated-Stirling matrix
  $A_{m,k}$, not on $O_n$.

### 2.2 Why this does not close (7)

(7) is a **root-velocity/sign** condition for $O_n(1,t)$ at a root $\mu$ of
$O_n$. The coefficients of $O_n$ are weighted sums over $m$ (13), with weights
$W(n,m,k)=n!(n-m+1)_k/m!$ depending on $n,m,k$. (H)+(17) constrain $A_{m,k}$
but do not by themselves imply:

- real-rootedness of $O_n(1,t)$ for all $n$ (a prerequisite for (7) as stated
  over real roots);
- interlacing $O_n \preceq O_{n+1}$;
- the residue sign $O_{n+1}(\mu)/O_n'(\mu)<0$;
- or the root-velocity bound (7).

A closure would require an additional implication of the form
$$
 (H),(17) \Longrightarrow \text{(7) or equivalent interlacing/Bezoutian positivity for }O_n.
$$
No such implication is proved here, and no mechanism in the existing notes
provides it. In particular, U-operator §7 explicitly disclaims closure of the
j=2 root-velocity or discriminant/resultant statements.

### 2.3 Why this does not close TP2/Bezoutian positivity for j=2

The Bezoutian positivity relevant to consecutive interlacing is
$M(O_n,O_{n+1})\succ0$ (equivalently all residues negative, equivalently (7)).
This is a polynomial-root/interlacing condition.

The TP2/ratio-monotonicity condition on falling-factorial coefficients is
$c_{n+1,i+1}c_{n,i}\ge c_{n+1,i}c_{n,i+1}$. Even if this were proved from
(H)+(17), §5.4 gives an exact obstruction: ratio monotonicity alone does not
force Bezoutian positive definiteness. The counterexample is:
$$
 P=10+11t+t^2,\quad Q=10+22t+100t^2,
$$
falling-factorial coefficient vectors $(10,12,1)$, $(10,122,100)$ with
strictly increasing ratios, but $\operatorname{Disc}_t(Q)=-3516<0$.
Therefore TP2/ratio monotonicity of $c_{n,i}$ is not sufficient for Bezoutian
positivity; an additional condition (ensuring real-rootedness / a stronger
Bezoutian certificate) is required.

Moreover, deriving $c_{n,i}$ TP2 from (H)+(17) is itself nontrivial: $c_{n,i}$
is obtained from $A_{m,k}$ by the two convolutions (13) then (18). The weight
$W(n,m,k)=n!(n-m+1)_k/m!$ depends on $n,m,k$ simultaneously, so the standard
2D TP2 Cauchy–Binet composition does not directly apply.

---

## 3. Exact remaining gap

The remaining gap after (H) is **not** (H) itself. It is one of the following
explicit targets (any one sufficient in its own route):

### Gap R1 (root-velocity / interlacing, direct)

Prove (7): for every root $\mu$ of $O_n(1,t)$,
$$
 (1+\mu)\mu_n'(1)>2\mu,
 \qquad\text{i.e.}\qquad
 \frac{O_{n+1}(1,\mu)}{O_n'(1,\mu)}<0.
$$
Equivalently prove (D)+(R) of j2-base §4 for all $x>0$.

### Gap R2 (falling-factorial Bezoutian, route closure)

Prove a condition on $c_{n,i}$ strong enough to imply $M(O_n,O_{n+1})\succ0$.
Ratio monotonicity $c_{n+1,i+1}c_{n,i}\ge c_{n+1,i}c_{n,i+1}$ is **insufficient**
by §5.4. The missing condition must therefore be stronger than 2-row TP2 of
$c_{n,i}$ (for example, a real-rootedness certificate for $O_n$ together with
ratio monotonicity, or a direct Bezoutian/TP positivity in the falling basis).

### Gap R3 (coefficient TP2 transfer, intermediate)

If attempting the transfer route: prove the 2-row TP2 minor
$$
 \Delta_{n,k}^{(d)}:=d_{n,k}d_{n+1,k+1}-d_{n,k+1}d_{n+1,k}\ge0
$$
for $d_{n,k}=[t^k]O_n(1,t)$, from (H),(17), and the explicit weight
$W(n,m,k)=n!(n-m+1)_k/m!$. Even this intermediate does not close R2 because
of the §5.4 obstruction; it would only close the ratio-monotonicity step.

---

## 4. Feasible falsification / finite evidence

Script: `check_route.py` (pure Python, exact `Fraction`).

Finite checks through $n\le 13$ ($n+1\le 14$):
- $c_{n+1,i}/c_{n,i}$ increasing in $i$: **all pass** (finite evidence).
- TP2 of $d_{n,k}=[t^k]O_n$: **all pass** (finite evidence).
- TP2 of $c_{n,i}$: **all pass** (finite evidence).

These are **finite evidence only** (AGENT_PROOF_POLICY §5): they do not prove
R1–R3. They also do not falsify the route; they indicate the remaining gap is
a proof gap, not an obvious small-$n$ sign failure.

A feasible falsification check for the "H+(17)=>R3" implication would be: treat
$A_{m,k}$ as nonnegative variables satisfying (H) and (17) and search for a
feasible assignment making some $\Delta_{n,k}^{(d)}<0$. This was **not**
performed here (it requires an explicit feasible-region search/optimization;
the actual $A_{m,k}$ are fixed). The fixed-$A_{m,k}$ finite check above found
no counterexample.

Direction/monotonicity audit for the §5.4 obstruction:
- $P=10+11t+t^2$, $Q=10+22t+100t^2$.
- Falling-factorial coeffs: $P=10(1)+12(t)_1+1(t)_2$, $Q=10(1)+122(t)_1+100(t)_2$.
- Ratios $122/12\approx10.166\ldots>1$, $100/1=100>10.166\ldots$: increasing. OK.
- $\operatorname{Disc}_t(Q)=22^2-4\cdot10\cdot100=484-4000=-3516<0$. OK.
- Conclusion: ratio monotonicity does not imply Bezoutian positive definiteness.
  Direction checked.

---

## 5. Conclusion

(H) is proved and removes one named bottleneck. However, (H)+(17) do **not**
close j=2 by themselves:

- They do not imply root-velocity (7): no proved implication from associated-Stirling
  coefficient inequalities to the $O_n$ root-velocity/residue sign.
- They do not imply TP2/Bezoutian positivity sufficient for j=2: even the
  natural falling-basis ratio-monotonicity target is insufficient for Bezoutian
  positive definiteness (5.4), and deriving that ratio monotonicity from
  (H)+(17) through (13),(18) is itself an unproved transfer.

The exact remaining inequality is R1, R2, or R3 above; R3 is weaker than R2,
and R2 is stronger than ratio monotonicity. No all-index closure is claimed.
