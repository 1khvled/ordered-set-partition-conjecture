# Transform preservation (T) for the original ordered big-block family

> **AUDIT STATUS — CONDITIONAL, NOT A PROOF.** An earlier audit banner
> incorrectly disputed the root-`0` statement. The correct values are
> `H_{0,j}(0)=1` and `H_{n,j}(0)=0` for every `n>=1`. The draft's root-`0`
> observation is therefore correct. Its Lemmas I–II remain conditional on the
> aligned all-positive-combination/common-interlacer property that Section 5
> calls open. Do not cite the conditional argument as an established proof.
> See `proof-final-transform-audit/NOTES.md` for the corrected audit.

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-final-transform` was created or edited.

## 1. Setting

Fix `j >= 2`. With

A_j(z) = sum_{r=1}^{j-1} z^r / r!,
C_j(z) = e^z - sum_{r=0}^{j-1} z^r / r!,
D_j(z) = 1/(1-A_j(z)),
g_j(z) = C_j(z)/(1-A_j(z)),

the decomposition of the structural track is

F_j(z,t) = 1/(1 - A_j - t C_j) = D_j(z)/(1 - t g_j(z)).

Writing

H_j(z,t) = 1/(1 - t g_j(z)) = sum_{r>=0} h_{r,j}(t) z^r,
H_{r,j}(t) = r! h_{r,j}(t),

and d_k = [z^k] D_j(z), the outer convolution is

O_{n,j}(t) = n! sum_{r=0}^n d_{n-r}^{(j)} h_{r,j}(t)
            = n! sum_{r=0}^n d_{n-r}^{(j)} H_{r,j}(t)/r!.

Equivalently, with the structural inner row H_{r,j}(t) = r! h_{r,j}(t),

O_{n,j}(t) = n! [z^n] (D_j(z) H_j(z,t)).                                    (T0)

Goal (T). Assuming (H) (each H_{r,j} is real-rooted with negative roots and
consecutive rows interlace), prove that for every n >= 0

- O_{n,j} is real-rooted, and
- O_{n,j} interlaces O_{n+1,j}.

## 2. Recurrence form of the transform

Since D_j = 1/(1-A_j), multiplying H_j by D_j gives

(1 - A_j(z)) sum_{n>=0} O_{n,j}(t) z^n / n! = sum_{r>=0} H_{r,j}(t) z^r / r!.

Extracting coefficients yields the exact recurrence (verified symbolically for
j=2,...,8, n<=12):

O_{n,j}(t) = H_{n,j}(t) + sum_{r=1}^{j-1} binom(n,r) O_{n-r,j}(t).          (T1)

(For n-r < 0 the term is zero.)

Thus (T) is equivalent to showing that recurrence (T1), whose forcing term is
the inner family H_{n,j} and whose feedback uses the previous j-1 rows,
preserves real-rootedness and interlacing.

## 3. Conditional proof via interlacing families

Let P_k(t) = O_{k,j}(t). Let q = floor(n/j). All polynomials P_{n-r}
(1 <= r <= j-1) and H_{n,j} have degree either q or q-1, because r < j.
(Verified for all tested j,n.)

Lemma I (interlacing family, window form). For every n >= j, the family

F_n = { H_{n,j} } union { binom(n,r) P_{n-r} : 1 <= r <= j-1 }

is an interlacing family: every member is real-rooted, and every positive
linear combination of the members is real-rooted.

Consequence of Lemma I (Marcus-Spielman-Srivastava). Taking the coefficients
1,...,1 in (T1) gives that P_n is real-rooted. This proves the first half of
(T) for all n, by induction on n, provided the induction hypothesis gives
real-rootedness of the P_{n-r} (true for r>=1 by induction) and of H_{n,j}
(true by (H)).

Lemma II (interlacing of consecutive rows). Under the same hypotheses,
P_n interlaces P_{n+1}.

Proof sketch of Lemma II. P_{n+1} is itself a positive combination of the
interlacing family

F_{n+1} = { H_{n+1,j} } union { binom(n+1,r) P_{n+1-r} : 1 <= r <= j-1 },

which contains P_n as the r=1 member. Since P_n interlaces every other
member of F_{n+1} (induction hypothesis: consecutive and window
interlacing), it interlaces any positive combination of them; in particular
it interlaces the sum of all members except itself. Adding back the
P_n-term with positive coefficient binom(n+1,1)=n+1 keeps the interlacing
by the standard Hermite-Biehler (Obreschkoff) lemma for p + lambda q when
p interlaces q. Hence P_n interlaces P_{n+1}.

Together, Lemmas I-II give (T) for all n, conditionally on Lemma I.

## 4. Computational evidence for Lemma I

For j=2,3,4,5, n<=12:

- Every member of F_n is real-rooted (follows from (H) and induction).
- All 200-300 random positive combinations of the members of F_n are
  real-rooted (sympy check). No counterexample found.
- Pairwise simple interlacing (alternating roots) holds for every pair in
  F_n for j=2,3,4,5, n<=12.
- Pairwise simple interlacing holds for the O-rows inside each window
  {P_{n-1},...,P_{n-j+1}} (0 failures for j=4, n<=12).

Scripts: test_interlace_family.py, test_interlace_family2.py,
test_interlace_family3.py, test_interlace_pairs.py, test_opairs.py,
test_recurrence2.py in this directory.

## 5. Precise gap

The proof above is conditional. What is NOT proved (and is the single
remaining obstruction for (T)) is:

Gap T1 (all-index interlacing family). Prove that for every j>=2 and every
n>=j, the family F_n in section 3 is an interlacing family. Equivalently,
prove that the roots of H_{n,j} and of each P_{n-r} interlace in a way that
is preserved under taking positive combinations, for the specific (non-PF)
kernels D_j.

This is not a finite check, and it does not follow from:

- Polya frequency / Branden-Leite (the kernel D_j is not PF for j>=3, and
  g_j is not PF; see proof-final-structural/NOTES.md section 2.1);
- planar-network / LGV (the natural network gives a reflected inequality,
  see proof-final-network/NOTES.md section 1);
- ordinary TN of the sheared Riordan array (blocked by Lemma 3 of
  proof-final-network/NOTES.md section 3).

What would close the gap: an induction proving that each P_{n-r} in the
window {n-j+1,...,n-1} interlaces H_{n,j} and that these interlacings are
aligned (same orientation) so that the whole family is interlacing. The
numerical data suggest the interlacing orientation is uniform (largest root
of the higher-degree members is 0 for H_{n,j}, all members have negative
roots otherwise), but a uniform-root-separation inequality analogous to the
weighted (WRS) lemma of proof-weighted/NOTES.md is still missing and appears
to be the exact quantitative content needed.

Status. (T) is proved conditionally on Gap T1. Gap T1 is open; all finite
evidence up to j<=5, n<=12 (and j<=8, n<=40 for real-rootedness alone) is
consistent with it.
