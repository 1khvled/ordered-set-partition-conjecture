# Lean 4 source

## Default build: verified and disk-safe

The default library target imports only the dependency-free core and builds
successfully:

```powershell
lake build
```

`OrderedBigBlocks/Core.lean` formalizes:

- binomial coefficients and integer coefficient-list polynomials;
- the first-block recurrence for `O_{n,j}`;
- the published initial `j=3` rows;
- the ordered-Bell specialization at `t=1`.

You can also check the core directly:

```powershell
lean OrderedBigBlocks/Core.lean
```

## Experimental Mathlib files

`OrderedBigBlocks/Basic.lean` and `OrderedBigBlocks/Sturm.lean` contain a more
ambitious Mathlib-based EGF and Sturm-certificate development. They have **not
compiled successfully** and are not imported by the default library target.
Mathlib is deliberately not a default dependency, so building the verified
core does not download or compile Mathlib.

The earlier local Mathlib cache was removed after it grew to 6.83 GB.
