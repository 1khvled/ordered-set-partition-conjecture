# Gap T-AUDIT attack: aligned common-interlacer propagation for the T1 window

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-common-interlacer` created/edited.  
Outcome: **no proof of Gap T-AUDIT obtained; no finite counterexample found.**  
This note gives a precise reduction of the gap to explicit boundary/common-interlacer
conditions, proves the easy implications, and records why the obvious closure routes fail.

## 1. Objects and checked identities

Fix `j>=2`. Let

- `A_j(z)=sum_{r=1}^{j-1} z^r/r!`, `C_j(z)=sum_{r>=j} z^r/r!`, `D_j(z)=1/(1-A_j(z))`.
- `g_j(z)=C_j(z)D_j(z)`, `H_j(z,t)=1/(1-t g_j(z))`, `F_j(z,t)=D_j(z)/(1-t g_j(z))`.
- `H_n(t)=n![z^n]H_j(z,t)`, `O_n(t)=n![z^n]F_j(z,t)`.

Proved identities (see `proof-final-transform-audit/NOTES.md`, reproved by its scripts):

1. `H_0=1`, `H_n(0)=0` and `t|H_n` for `n>=1`, `deg_t H_n=floor(n/j)`.
2. `O_n(0)=n![z^n]D_j>0`, `deg_t O_n=floor(n/j)`, nonnegative coefficients.
3. **(T1)** `O_n = H_n + \sum_{r=1}^{j-1} \binom{n}{r} O_{n-r}`, `n>=0`.

Define the window

`W_n = {H_n} \cup {B_{n,r}(t) := \binom{n}{r} O_{n-r}(t) : 1<=r<=j-1},  n>=j`.

Degree split with `q=floor(n/j)`, `s=n mod j`:
`deg H_n=q`, `deg B_{n,r}=q` for `1<=r<=s`, `deg B_{n,r}=q-1` for `s<r<=j-1`.

## 2. Common interlacer: exact criterion

For real-rooted polynomials, use the standard mixed-degree interlacing convention:
`f` interlaces `g` if `deg g=deg f+1` and roots alternate
`g1<=f1<=g2<=...<=f_d<=g_{d+1}`. Write `f << g` for this relation.
A family `F` of real-rooted polynomials has a **common interlacer** `G` if every
nonzero member of `F` interlaces `G`.

### Lemma 2.1 (common interlacer => real-rooted positive combinations)

Let `F={f_1,...,f_k}` be real-rooted polynomials of degrees in `{q,q-1}` with
positive leading coefficients, and let `G`, `deg G=q+1`, be a common interlacer.
Then every positive combination `P=sum_i lambda_i f_i`, `lambda_i>0`, is
real-rooted. In particular (T1) gives `O_n` real-rooted if `W_n` has a common
interlacer and all its members are real-rooted.

*Proof sketch.* This is the standard Dedieu/Chudnovsky-Seymour/Fell implication
used in MSS interlacing families. For each root interval of `G`, the sign of a
degree-q member at the endpoints is forced by interlacing; a positive sum has the
same forced sign pattern and therefore has a root in every interval. Degree-q-1
members are treated by the same interval argument after viewing them through the
common degree-q+1 interlacer. A complete proof requires the usual sign/degree
normalization; this implication is the easy direction.  ∎

### Lemma 2.2 (pairwise interlacing is insufficient)

For families of size >=3, pairwise interlacing does **not** imply a common
interlacer. Therefore finite pairwise checks cannot prove (CI_n).

*Reason.* The classical obstruction is a 3-cycle in the root-order compatibility
graph (Fell/Dedieu). Equivalently, the cone of positive combinations may leave the
real-rooted locus even when every pair interlaces.

## 3. The propagation target

The audit asks for an induction step `W_n -> W_{n+1}`. The next window is

`W_{n+1} = {H_{n+1}} \cup {B_{n+1,1}O_n} \cup {B_{n+1,r}O_{n+1-r} : 2<=r<=j-1}`.

The carried-over tail is `T_{n+1}={O_{n-1},...,O_{n-j+2}}` (up to positive
scaling). Thus propagation requires four simultaneous alignments:

- **(P1) Common interlacer for carried tail.** Needed as part of (CI_{n+1}).
- **(P2) O_n aligns with carried tail.** `O_n` must interlace
  `O_{n-1},...,O_{n-j+2}` in the same orientation/common-interlacer frame.
- **(P3) H_{n+1} aligns with O_n and carried tail.** Cross-family boundary:
  the new inner polynomial must join the O-window.
- **(P4) Degree-shift coherence.** When `n+1 ≡ 0 (mod j)`, degrees in the
  window move from `q,q-1` to `q+1,q`; the common interlacer degree becomes
  `q+2`.

### Proposition 3.1 (exact reduction)

Assume (H): every `H_m` is real-rooted and `t|H_m` for `m>=1`.
Assume also (O-chain): `O_0,...,O_n` are real-rooted and consecutive
`O_{k-1},O_k` interlace in the standard sense for all `k<=n`.
Then `W_{n+1}` has a common interlacer **if and only if** the two boundary members
`O_n` and `H_{n+1}` can be placed into a common interlacer of the O-tail:

> there exists `G_{n+1}` with `deg G_{n+1}=floor((n+1)/j)+1` such that
> every element of `{H_{n+1}, O_n, O_{n-1}, ..., O_{n-j+2}}` interlaces
> `G_{n+1}` after the appropriate degree normalization.

*Proof.* `W_{n+1}` is exactly this set up to positive scalar multiples, which do
not affect interlacing. Conversely any common interlacer of `W_{n+1}` is a common
interlacer of this normalized set. The O-chain assumption supplies pairwise
compatibility along the tail but, by Lemma 2.2, not a common interlacer for the
whole enlarged family.  ∎

This isolates the gap: the O-chain alone is insufficient; the two new boundary
members must be simultaneously compatible.

## 4. Why the obvious routes do not close the gap

1. **Pairwise interlacing induction.** Proving every adjacent pair interlaces
   gives P1/P2 pairwise data but not a common interlacer for size>=3. The T1
   sum then is not forced real-rooted by pairwise data alone.
2. **Convex-combination / MSS random-combination evidence.** Finite random
   positive combinations (including `diag_window.py`) are evidence only; they do
   not certify the whole cone.
3. **Wronskian sign.** For two polynomials, constant-sign Wronskian is the right
   certificate. For a whole window, one needs a common interlacer/Wronskian matrix
   condition (e.g. Hermite-Biehler/Bezoutian or a common G) simultaneously for
   all members; pairwise Wronskian signs do not compose.
4. **H/O cancellation at 0.** Since `H_n(0)=0` and `O_n(0)>0`, the boundary at 0
   is an orientation constraint, not a cancellation obstruction. It forces the
   top root of O_n to lie strictly below 0 in any compatible interlacing.
5. **j=2 reduction.** For j=2, `W_n={H_n, nO_{n-1}}` has size 2, so common
   interlacing is pairwise. The propagation still needs `H_{n+1}` to interlace
   `O_n`, a cross-family boundary step that is not proved here.
6. **Non-PF kernels.** `g_j` and `D_j` are not PF in the relevant ranges
   (`proof-final-structural`, `proof-final-transform-audit`), so PF/TN multiplier
   theorems do not supply a common interlacer.

## 5. Finite diagnostics (evidence only)

Scripts in this directory:

- `diag_window.py` - pairwise interlacing, random positive-combination real-rootedness,
  and boundary alignment checks for `W_n -> W_{n+1}`.
- `diag_wronskian.py` - Wronskian real-root/sign checks between `H_n` and window members.
- `diag_search.py` - larger random search.

Observed (not a proof): for `j=2..6`, `n<=19` (as available in `diag_search.out`/`diag_window.out`),
no pairwise failure, no random positive-combination failure, and no propagation
pairwise failure was found. Wronskians tested had no detected real roots on the
sampled range. This supports the conjectured propagation but does not close Gap T-AUDIT.

## 6. Exact obstruction

### 6.1 Geometric common-interlacer condition

Normalize every nonzero member of `W_n` to positive leading coefficient. Let
`q=floor(n/j)`. A degree-`q+1` common interlacer `G_n` with simple roots
`g_1<...<g_{q+1}` exists iff the following interval constraints are simultaneously
satisfiable:

- For a degree-`q` member `f` with roots `f_1<=...<=f_q`, interlacing requires
  `g_i <= f_i <= g_{i+1}`.
- For a degree-`q-1` member `f`, interlacing requires `g_i <= f_i <= g_{i+2}`
  in the standard mixed-degree embedding.
- For `H_n` (degree `q`, `H_n(0)=0`, `t|H_n`), the largest root is `0`, forcing
  `g_q <= 0 <= g_{q+1}`.
- For every `O_k` member, `O_k(0)>0` forces its largest root to be strictly
  negative; compatibility with `H_n` then forces that largest root to lie in
  `(g_q,0)`.

Thus the gap is an all-index nonemptiness of this intersection-of-root-interval
system, together with propagation of that nonemptiness from `n` to `n+1`.

### 6.2 j=2 boundary obstruction

For `j=2`, `W_n={H_n,n O_{n-1}}` has size 2, so (CI_n) is equivalent to pairwise
interlacing of `H_n` and `O_{n-1}`. The propagation step is still open:

> prove `H_{n+1}` interlaces `O_n = H_n + n O_{n-1}` for all `n>=2`.

This is a concrete cross-family boundary inequality; it is weaker than the full
common-interlacer propagation for `j>=3` but not implied by the checked
identities.

### 6.3 All-index statement

> **Gap T-AUDIT (aligned common interlacer propagation).** Fix `j>=2`.
> Assume (H) and the T1 identities. Prove or disprove: for every `n>=j`, the
> window `W_n={H_n} \cup {binom(n,r)O_{n-r}:1<=r<=j-1}` has a common interlacer
> `G_n`, and `G_n` can be chosen so that `O_n=H_n+sum_r binom(n,r)O_{n-r}` is
> aligned as a carried member of `W_{n+1}` together with `H_{n+1}`.

Equivalent missing certificates: (i) a simultaneous common interlacer for
`{H_{n+1}, O_n, ..., O_{n-j+2}}`; (ii) a linear-closure proof that every positive
combination in `W_n` is real-rooted, propagated through degree shifts; or (iii)
another all-index transform theorem converting (H) to (T) for this non-PF kernel.
None is established here.

## 7. Files

- `NOTES.md` - this report.
- `diag_window.py`, `diag_window.out` - window/pairwise/propagation diagnostics.
- `diag_wronskian.py`, `diag_wronskian.out` - Wronskian diagnostics.
- `diag_search.py`, `diag_search.out` - random search.
