# Lean 4 source

This directory contains a Mathlib-based formalization of:

- the first-block recurrence for `O_{n,j}`;
- the exact EGF equation `(1 - B_j) F_j = 1`;
- an executable rational Sturm chain;
- a finite real-root certificate based on disjoint rational sign-change intervals.

## Build status

`OrderedBigBlocks/Core.lean` is dependency-free and **compiles successfully**
with the installed Lean 4 toolchain:

```powershell
lean OrderedBigBlocks/Core.lean
```

It formalizes the coefficient recurrence and verifies the initial `j = 3` rows
and the ordered-Bell specialization.

The Mathlib files `Basic.lean` and `Sturm.lean` **have not compiled
successfully**. The initial build attempt reported elaboration errors, and the
local Mathlib dependency/build cache was removed after it grew to 6.83 GB.
Those files are included as a transparent starting point, not as a verified
formal proof of the conjecture.

A future build will download and compile Mathlib again and may require several
gigabytes of free disk space:

```powershell
cd lean
lake build
```

Do not run this on a low-disk machine unless that space is available. The pinned
dependency revisions remain in `lake-manifest.json`.
