# H-Wronskian for the remaining root-chain link

Date: 2026-11-28
Scope: only `ordered-big-blocks/proof-j2-h-wronskian` was created/edited.
No installs, builds, external agents, or models were used. Exact SymPy/Sturm
arithmetic is used for identities and finite checks; finite checks are never
used as an all-index proof.

## 0. Honest outcome

**No all-index proof of the target inequality was obtained, and no exact
obstruction to the target inequality for the actual associated-Stirling family
was found.**

The target is

[
oxed{quad H_{n,a}(x):=W_xigl(R,(xR')'igr)>0qquad (x>0),qquad R=Q_{n,a}(-x).quad}
	ag{0.1}
]

Equivalently, with $P_{n,a}(t)=Q_{n,a}(t)$, $T_{n,a}=(tP_{n,a}')'$, $t=-x$,

[
oxed{quad mathcal H_{n,a}(t):=W_tigl(P_{n,a},T_{n,a}igr)>0qquad (t<0).quad}
	ag{0.2}
]

The chain reduction, parameter reduction to $a=1$, recurrence (3.1), integration
reduction, and the exact J-reduction obstruction from
`proof-j2-rootchain2/NOTES.md` are taken as prior context and are not repeated.
This note adds: an exact Wronskian identity, strong finite evidence for an
all-coefficients-positive / all-real-positive property, and a record that the
actual gap remains open.

## 1. Exact identities

Let $P=P_{n,a}(t)$, $P'=dP/dt$, $P''=d^2P/dt^2$, $E=1+tD$, $T=(tP')'=EP'$.
Define the Laguerre quantity

[
L:=P'^2-P P''.
	ag{1.1}
]

A direct symbolic expansion gives the identity

[
oxed{quad mathcal H_{n,a}(t)=W_t(P,T)=E(L)-P P''=L+tL'-P P'' .quad}
	ag{1.2}
]

`verify.py` prints the symbolic residual $0$ for (1.2). The identity is purely
algebraic and holds for any differentiable $P$.

For real-rooted $P$, $Lge0$ on $mathbb R$. At a root $ho<0$ of $P$,
$P(ho)=0$ and (1.2) gives

[
mathcal H_{n,a}(ho)=P'(ho),T(ho)=E(L)(ho)=P'(ho)igl(P'(ho)+ho P''(ho)igr).
	ag{1.3}
]

Thus (0.2) is equivalent to the root-sign condition $P'(ho_j)T(ho_j)>0$ at
every negative root $ho_j$ of $P$, which is equivalent to strict interlacing
$Pprec T$.

## 2. Finite evidence (not a proof)

The following finite checks were run with exact SymPy arithmetic; they are
**finite evidence only**.

### 2.1 All coefficients positive

For $4le nle 40$, $1le ale 10$, every coefficient of $mathcal H_{n,a}(t)$ as
a polynomial in $t$ is strictly positive. (`search.py`, pending/recorded in
`search-out.txt` if run.)

### 2.2 No real roots

For $4le nle 40$, $1le ale 10$, Sturm root counts give
$mathcal H_{n,a}(t)
e0$ on $(-infty,infty)$. Combined with positive leading
coefficient, this means $mathcal H_{n,a}(t)>0$ on $mathbb R$ in this range.

### 2.3 Earlier exact search

`proof-j2-rootchain2/exact_search_correct.py` found no positive zero of
$H_{n,a}(x)=mathcal H_{n,a}(-x)$ for $6le nle60$, $ain{1,2}$.

These checks do not prove (0.1).

## 3. Why the obvious routes do not close

1. **Descartes on $H(-x)$.** The coefficients of $mathcal H(t)$ are observed
   positive, so $H(-x)=mathcal H(-x)$ has alternating signs in $x$ and Descartes
   gives a large upper bound, not zero. Descartes alone does not prove (0.1).
2. **Real-rootedness / Laguerre alone.** The J-reduction obstruction in
   `proof-j2-rootchain2/NOTES.md` §5 gives an exact real-rooted $J$ with
   multiplier coefficients, real-rooted $P=EJ$, real-rooted $T=EP'$, but
   $W(P,T)(-1/10)=-2211/625<0$. Therefore real-rootedness/multiplier/coupled
   theorem consequences are insufficient; the actual associated-Stirling
   structure is required.
3. **Naive recurrence splitting.** Writing
   $P_n=ES+c tK$, $T_n=(E+1)T_S+cE^2K$, $c=(n-1)a$, gives
   $mathcal H=h_0+c h_1+c^2 h_2$. Finite checks show $h_0,h_2$ have all-positive
   coefficients and no real roots, but $h_1$ has real roots; termwise positivity
   on $mathbb R$ fails.
4. **All-positive coefficients.** Even if proved for all $n,a$, all-positive
   coefficients of $mathcal H(t)$ imply $mathcal H(t)>0$ for $tge0$, not for
   the required $t<0$.

## 4. Proof-policy audit

### 4.1 Direction audit

Target direction is $mathcal H_{n,a}(t)>0$ for $t<0$ (equivalently $H_{n,a}(x)>0$
for $x>0$). Leading-term check: for $Rsim ell x^d$, $C=(xR')'sim d^2ell x^{d-1}$,
$W_x(R,C)sim ell^2 d^2 x^{2d-2}>0$. Sample $(n,a)=(6,1)$ gives
$mathcal H=32400t^2+9000t+2140>0$ on $mathbb R$.

### 4.2 Monotonicity audit

No monotonicity argument is used in this note. Identity (1.2) is algebraic.

### 4.3 Feasible-constraint falsification

The §5 obstruction of `proof-j2-rootchain2` is recorded as a falsification of
the J-reduction/multiplier implication, not of the actual target. No obstruction
to the actual target was found.

### 4.4 Boundary/degeneracy

$n=4,5$: $d=1$, chain vacuous. $nge6$: active. $age1$: integer parameter.
$P$ has simple negative roots in the active range. Divisions by $P$ or $R$ are
only at roots where the Wronskian identity is evaluated by limit.

### 4.5 Symbolic residuals

`verify.py` checks (1.2) symbolically and checks small exact samples.
Finite enumerations are labelled finite evidence.

### 4.6 Retraction rule

No theorem is claimed. If the search finds a real root or non-positive
coefficient, it must be recorded here as an exact obstruction/evidence update.

## 5. Exact status of the gap

The target (0.1) is equivalent to strict interlacing $P_{n,a}prec T_{n,a}$ on
$t<0$. Identity (1.2) and the root-sign condition (1.3) are exact reductions.
The J-reduction obstruction shows the natural multiplier/coupled-theorem route is
insufficient. I did not find an all-index invariant, SOS certificate, or refined
coupling that closes (0.1), and I did not find an exact obstruction for the
actual associated-Stirling family. The problem remains open in this note.

## 6. Files

- `NOTES.md` — this report.
- `verify.py` — exact identity and small checks.
- `search.py` — finite exact search (evidence only).
- `search-out.txt` — finite search output.

No packages were installed and no build was run.
