# Ordered set partitions counted by big blocks

![Exact regression tests](https://github.com/1khvled/ordered-set-partition-conjecture/actions/workflows/tests.yml/badge.svg)

A reproducible exact-arithmetic attack on Per Alexandersson's 2026 conjecture

\[
O_{n,j}(t)=\sum_{\pi}t^{\#\{B\in\pi:|B|\ge j\}},\qquad
\sum_{n\ge0}O_{n,j}(t)\frac{z^n}{n!}=\frac1{1-B_j(z,t)},
\]

where the sum is over **ordered** set partitions and

\[
B_j(z,t)=\sum_{r=1}^{j-1}\frac{z^r}{r!}
+t\left(e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}\right).
\]

The conjecture asserts that every \(O_{n,j}\) is real-rooted and that
\(O_{n,j}\interl O_{n+1,j}\).

> **Status:** research in progress. The universal differential recurrence is
> proved, and the all-parameter real-rootedness/interlacing statement remains
> open. The exact remaining consecutive-interlacing target is
> \(O_{n+1,j}(\mu)/O_{n,j}'(\mu)<0\) at every root \(\mu\) of \(O_{n,j}\);
> the lag-\(j\) WRS route is equivalently a negative-axis Wronskian inequality.
> Finite exact checks are evidence, not a proof of the full conjecture.
> See [RESULTS.md](RESULTS.md) for the concise proved/computed status.

## Main proved recurrence and proof obstruction

The universal recurrence is
\[
O_{n,j}=nO_{n-1,j}+2T O_{n-1,j}
+\sum_{r=2}^{j-1}\binom nr O_{n-r,j}
+\binom{n-1}{j-1}t[1+(t-1)\partial_t]O_{n-j,j},
\quad T=t\partial_t.
\]
A naive compatible-window induction fails: compatibility of the lagged row
with the main part is not enough, because the derivative term does not vanish
at the lagged row's zeros. A rigorous weighted root-separation inequality is
given in `proof-weighted/NOTES.md`, together with an explicit obstruction to
the unweighted version.

## Correct coefficient recurrence

Choosing the first block gives the integer recurrence

\[
O_{0,j}=1,
\qquad
O_{n,j}=
\sum_{r=1}^{j-1}\binom nr O_{n-r,j}
+t\sum_{r=j}^{n}\binom nr O_{n-r,j}.
\]

The binomial coefficient is \(\binom nr\), **not** \(\binom{n-1}{r-1}\):
the latter forces the block containing \(1\) to be first and undercounts
ordered partitions.

For example, for \(j=3\),

```text
O_0 = 1
O_1 = 1
O_2 = 3
O_3 = 12 + t
O_4 = 66 + 9t
O_5 = 450 + 91t
O_6 = 3690 + 973t + 20t^2
```

## Repository layout

- `scan.py`, `interlace.py`, `frontier.py` — exact coefficient generation,
  certified root isolation, Sturm cross-checks, discriminants and Newton tests.
- `guess_recurrence.py` — exact search for short differential recurrences.
- `numerical/` — independent corrected numerical audit and machine-readable logs.
- `proof/` — proof-track identities, candidate recurrences and proof gaps.
- `proof-weighted/` — rigorous weighted root-separation theorem and exact obstruction.
- `proof-matrix/` — Riordan/shear factorization and an exact negative-minor obstruction.
- `proof-bezout/` — Bézoutian decomposition and exact indefinite-kernel obstruction.
- `proof-final-*` — final targeted attacks: corrected WRS/Wronskian, consecutive
  sign conditions, inner PDE, threshold renewal, middle zeros, crossing
  involutions, and common-interlacer propagation.
- `proof-global-wronskian/`, `proof-all-lag-pick/`, `proof-cross-boundary/`,
  `proof-residue-domination/`, `proof-joint-interlacing/` — exact reductions
  showing that the Pick, boundary, and joint-interlacing routes are equivalent
  to or harder than the original problem; no solution is claimed.
- `proof-chapter-renewal/`, `proof-branching-coupling/` — positive Möbius flow
  and branching-coupling constructions, with the exact non-positive
  coefficient-extraction obstruction.
- `proof-j2-base/`, `proof-j2-resultants/`, `proof-associated-stirling/`,
  `proof-associated-stirling-lit/`, `proof-associated-stirling-slc/`,
  `proof-factorial-row/` — the \(j=2\) deformation, resultant equivalence,
  literature check, and the associated-Stirling reduction of the remaining
  ratio inequality to a single cross term or to one operator-preservation lemma.
- `proof-j2-highn/` — new certified high-row search for \(j=2\), disjoint from
  the other counterexample hunts.
- `proof-automated-consecutive/`, `proof-hurwitz-wronskian/` — exact finite
  discovery of the consecutive-Wronskian invariant and its Routh/Hermite
  obstructions; numerical left-half-plane evidence is not used as a certificate.
- `proof-counterexample-new/`, `proof-large-j-sweep/` — new-range counterexample
  hunts with compact exact logs. Large-threshold results are reported from the
  checkpoint only until the full sweep finishes.
- `proof-shift/` — exact compatibility probes and the audit trail for a retracted lemma.
- `literature.md` — sourced context and distinction from related ordered lifts.
- `lean/` — a compiling dependency-free Lean recurrence plus Mathlib-based
  Sturm-certificate source that is not yet verified.
- `paper/` — manuscript source and build instructions.

## Python environment

Tested with CPython 3.14.4, `python-flint` 0.9.0 and `gmpy2` 2.3.2.

```powershell
py -3.14 -m pip install python-flint gmpy2 sympy
py -3.14 scan.py --jmin 3 --jmax 20 --nmax 150 --out scan-results.jsonl
py -3.14 interlace.py --jmin 3 --jmax 20 --nmax 150 --wronskian
py -3.14 frontier.py --small-nmin 141 --small-nmax 600 --jmin 6 --jmax 60
py -3.14 numerical/compact_frontier.py --mode wide --prec 192 --out numerical/compact-wide.jsonl --summary numerical/compact-wide-summary.json
```

All final root decisions use rigorous Arb enclosures or exact Sturm arithmetic;
float64 root finders are not used as certificates.

## Lean 4

The default Lean build is dependency-free and verified:

```powershell
cd lean
lake build
```

Experimental Mathlib EGF/Sturm files are present but unverified and excluded
from the default target.

## References

- Per Alexandersson, *Open problems on real-rootedness*, SymCat (2026):
  <https://www.symmetricfunctions.com/realRootedOpenProblems.htm#conjOrderedBigBlocks>
- Public computational audit, Issue #108:
  <https://github.com/1337isnot1337/solved-conjectures/issues/108>
