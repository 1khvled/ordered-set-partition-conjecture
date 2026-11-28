# j=2 B-split interlacing target: rigorous reduction and exact obstructions

Date: 2026-11-28. Scope: only `proof-j2-b-split` created/edited. No installs,
no builds, no other agents/models.

## Outcome

**No unconditional all-`n` proof of the target interlacing was obtained.**

The target from `proof-new-split/NOTES.md` §5 holds in all exact finite checks
performed here (`n<=25`), but proving it for all `n` would close the open `j=2`
consecutive-interlacing / real-rootedness conjecture (SymCat issue #108; see
`proof-j2-resultants/NOTES.md` Lemma 1–2 for the equivalent x-deformation
obstruction). The natural closure routes are blocked by the exact obstructions
recorded below. Finite checks are evidence only.

---

## 1. Setup and exact identities (all confirmed)

Let `j=2`, `Q_n := O_{n,2}`, `Q_0=1`, `deg Q_n = \lfloor n/2\rfloor`,
`Q_n(0)=n!`, all checked roots real, simple, negative. Define

```
A_n = n Q_{n-1} + (n-1) t Q_{n-2},
B_n = 2 Q'_{n-1} + (n-1)(t-1) Q'_{n-2}.                     (1)
```

Then `Q_n = A_n + t B_n`. The derivative-free identity is exact:

```
B_n = \sum_{r=2}^{n-1} \binom{n-1}{r} (Q_{n-r} + Q_{n-1-r}),  n>=2.  (2)
```

Hence `B_n` has nonnegative coefficients for `n>=3` (a positive combination of
nonnegative-coefficient polynomials). Verified exactly for `n<=25`.

---

## 2. Equivalent forms of the target

Target: for all `n>=6`, `Q_{n-1}` interlaces `t B_n`.

Let `d=\deg Q_{n-1}=\lfloor(n-1)/2\rfloor`. Then `\deg B_n=d-1`, `\deg(tB_n)=d`.
The Wronskian `W(f,g)=f'g-fg'` satisfies

```
W(Q_{n-1}, tB_n)(mu) = mu * Q'_{n-1}(mu) * B_n(mu)            (3)
```

at every root `mu<0` of `Q_{n-1}`. Since `mu<0` and `Q'_{n-1}(mu)!=0`, the
no-real-zero/constant-sign condition `W(Q_{n-1},tB_n)<0` is equivalent to

```
Q'_{n-1}(mu) * B_n(mu) > 0   at every root mu<0 of Q_{n-1}.   (4)
```

Substituting (1) gives `B_n(mu)=2Q'_{n-1}(mu)+(n-1)(mu-1)Q'_{n-2}(mu)`, so (4)
is algebraically equivalent to the ratio form

```
R_n(mu) := (n-1)(1-mu) * Q'_{n-2}(mu) / Q'_{n-1}(mu) < 2.     (5)
```

Also equivalent is `B_n` interlacing `Q_{n-1}` as a degree `(d-1,d)` pair:
`W(B_n,Q_{n-1})(mu) = -Q'_{n-1}(mu)*B_n(mu)`, hence
`W(Q_{n-1},tB_n)(mu) = -mu * W(B_n,Q_{n-1})(mu)` and the two Wronskians have the
same sign at roots of `Q_{n-1}`. Both have constant negative sign at 0 in all
checks (`n<=25`).

---

## 3. Rigorous cone reduction (assuming induction)

Assume the induction hypothesis `Q_{n-2} <= Q_{n-1}` (consecutive interlacing,
all roots real/simple/negative). Then:

- `Q_{n-1}` interlaces `n Q_{n-1}` trivially.
- `Q_{n-1}` interlaces `t Q_{n-2}`: at a root `mu` of `Q_{n-1}`,
  `W(Q_{n-1},tQ_{n-2})(mu)=mu*Q'_{n-1}(mu)*Q_{n-2}(mu)`. By interlacing of
  `Q_{n-2},Q_{n-1}` with positive leading coefficients,
  `Q_{n-2}(mu)/Q'_{n-1}(mu)>0`, and `mu<0`, so `W<0`. Also
  `W(Q_{n-1},tQ_{n-2})(0)= -Q_{n-1}(0)Q_{n-2}(0)= -(n-1)!(n-2)!<0`.
- By Obreschkoff's cone, `Q_{n-1}` interlaces `A_n = nQ_{n-1}+(n-1)tQ_{n-2}`.

Therefore, given `Q_{n-2}<=Q_{n-1}`, the target (4)/(5) implies
`Q_{n-1} <= Q_n = A_n+tB_n`. This is a rigorous sufficient reduction. The base
cases `n<=5` are exact-checked.

---

## 4. Finite evidence (exact, not a proof)

`verify.py` checks, for `n<=25`:

1. B identity (1)=(2) and `B_n>=0` coefficients.
2. `Q_{n-1}` interlaces `tB_n` (`W` has no real zeros, `W(0)<0`).
3. `B_n` interlaces `Q_{n-1}` (equivalent form).
4. `R_n(mu)<2` at every negative root of `Q_{n-1}`.
5. `Q_{n-1}` interlaces `A_n`.
6. Individual B-sum terms `C(n-1,r)(Q_{n-r}+Q_{n-1-r})(mu)/Q'_{n-1}(mu)` can be
   negative; cancellation is required.
7. `R_n` at the largest root (closest to 0) is increasing on the checked range
   (finite evidence only): `n=4: 0.795918`, ..., `n=25: 1.951517`.
8. `B_n` interlaces `B_{n-1}`.

All pass.

---

## 5. Exact obstructions

### 5.1 Termwise positivity fails

The B-sum (2) is coefficient-nonnegative but not termwise positive at roots of
`Q_{n-1}`. For example `n=6`, largest-magnitude root has two negative summand
ratios out of four; the total sum remains positive. Therefore no proof can
deduce (4) by showing every summand of (2) has the same sign as `Q'_{n-1}(mu)`.

### 5.2 Sharpness / vanishing margin (finite evidence)

On the checked range `n<=25`, `R_n(mu*)` at the largest root `mu*` of `Q_{n-1}`
increases toward 2 (`1.951517` at `n=25`). This is finite evidence only; no
all-`n` monotonicity or limit is proved. It indicates vanishing margin on the
checked range, so any all-`n` proof would need to be sharp at the hardest root
and cannot rely on uniform slack alone.

### 5.3 Equivalence in strength to an open conjecture

By §3, an all-`n` proof of (4)/(5) would prove `Q_{n-1}<=Q_n` for all `n>=2`,
i.e. the `j=2` consecutive real-rootedness/interlacing conjecture. This is open
(SymCat issue #108; `proof-j2-resultants/NOTES.md` proves the x-deformed
`(D)+(R)` statements are equivalent to this open problem). The present B-split
is a sufficient condition, not shown equivalent to the full x-deformation, but
closing it has the same open-conjecture consequence.

### 5.4 No closed residue/cancellation identity found

The ratio form (5) requires controlling `Q'_{n-2}/Q'_{n-1}` at roots of
`Q_{n-1}`. Derivatives `Q'_{n-2},Q'_{n-1}` interlace (Rolle, checked `n>=5`),
but at roots of `Q'_{n-1}`, not `Q_{n-1}`. The B-sum tail reappears at roots of
`Q_{n-1}` after eliminating `Q'_{n-1}` through the recurrence, and no
sign-closed residue identity or two-step invariant was found that makes the
required cancellation manifest for all `n`.

---

## 6. Policy audit (AGENT_PROOF_POLICY.md)

No proof is claimed. The following audits were completed (`audit.py`):

1. **Direction audit.** Sample `n=6, mu=-0.461030` (largest root of `Q_5`):
   `Q'_5(mu)=209.573376`, `B_6(mu)=130.829710`, product `27418.424053 > 0`;
   `R_6(mu)=1.375733 < 2`. The identities `W(Q_{n-1},tB_n)(mu)=mu Q'_{n-1}(mu)B_n(mu)`
   and `W(B_n,Q_{n-1})(mu)=-B_n(mu)Q'_{n-1}(mu)` were checked as symbolic
   residuals at roots for `n=4..9`. Cone-reduction signs: `Q_{n-2}(mu)/Q'_{n-1}(mu)>0`
   and `W(Q_{n-1},tQ_{n-2})(mu)<0` checked at all negative roots, `n=4..11`;
   `W(Q_{n-1},tQ_{n-2})(0)=-(n-1)!(n-2)!` checked `n=3..11`.
2. **Monotonicity audit.** No all-`n` monotonicity is claimed. The observation
   "`R_n(mu*)` increasing on checked range `n=4..25`" is finite evidence only.
3. **Feasible-constraint falsification.** Coefficient-nonnegativity of `B_n`
   alone is insufficient: witness `mu=-2.366243` (root of `Q_5`),
   `Q'_5(mu)=-209.573376`, `B_tilde(mu)=5600.104132` for
   `B_tilde=1+1000*t^2` (same degree, nonnegative coeffs), product
   `-1173632.729813 < 0`. Therefore any proof must use the specific B-split
   structure/cancellation, not just `B_n>=0`.
4. **Boundary/degeneracy.** Checked `n=0..7`: `deg Q_n=floor(n/2)`,
   `Q_n(0)=n!`, positive leading coefficient; degree formula checked through
   `n=15`. The small-`n` base cases `n<=5` are exact-checked in `verify.py`.
5. **Symbolic residuals.** `B_direct-B_sum==0` and `Q_n-(A_n+tB_n)==0` verified
   for `n<=15` in `audit.py`; also `n<=25` in `verify.py`. All finite
   enumerations are labeled finite evidence.

## 7. Files

- `NOTES.md` — this report.
- `verify.py` — exact finite checks (SymPy/Sturm; finite evidence only).
- `audit.py` — policy audit (directions, residuals, boundaries, falsification).
- `explore.py` — small exploratory script used during the search.

Run from this directory:

```powershell
py -3.14 verify.py
py -3.14 audit.py
```
