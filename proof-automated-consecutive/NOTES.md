# Automated consecutive sign discovery

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-automated-consecutive`.

## Status

No all-index proof or exact obstruction was found by the automated search. The
search strongly supports the consecutive sign condition through a global
Wronskian positivity phenomenon.

Verified findings (exact Sturm / high-precision evaluation):

1. For `j=2, n<=34`; `j=3, n<=24`; `j=4, n<=19`; `j=5, n<=19`,
   `W(O_n,O_{n+1})(t)` has no real roots and positive leading coefficient.
   Equivalently it is strictly positive on `R`.
2. For the same ranges, evaluating `Rho_{n+1}` at every root of `O_n` gives a
   strictly negative value (worst margin not attained at 0).
3. Every tested `W(O_n,O_{n+1})` has all coefficients positive. High-precision
   root enclosures suggest that its complex roots have negative real parts
   (Hurwitz stability) for the checked cases; this numerical observation is not
   used as an exact Hurwitz certificate.

Consequence for checked ranges: by induction from `O_0=1`, `O_n` is real-rooted
and `O_{n-1} interl O_n` for all checked `n`.

## Discovered candidate invariant

Candidate invariant (not proved):

> For all `j>=2, n>=0`, `W_n(t)=O_{n,j}(t)O_{n+1,j}'(t)-O_{n,j}'(t)O_{n+1,j}(t)`
> is a Hurwitz polynomial: all its roots lie in the open left half-plane.

This is stronger than the required consecutive sign condition and is consistent
with all exact checks. If proved, it implies real-rootedness and consecutive
interlacing by induction.

## What was tried

- Sum-of-squares decomposition search: `W_n` is positive, but no uniform
  SOS/Bezoutian factorization with simple coefficient recurrences was found.
- Coefficient recurrence search: no short recurrence with polynomial-in-`n`
  coefficients for the `W_n` coefficient sequences was found in the tested
  ranges.
- Residue sign search: no positive `Rho` value was found; no exact
  counterexample to domination was found.

## Remaining gap

Prove the Hurwitz stability of `W_n` for all `n,j`, or at least the weaker
global positivity of `W_n` on the real axis. The recurrence for `O_n` mixes Christoffel transforms of recent
rows with a lagged `tQ` term; the latter is not individually positive, so a
uniform proof requires a global cancellation not yet established.

## Files

- `automated_consecutive.py` - exact Sturm certification and Rho evaluation for
  a chosen range.
- `Wsturm2.py`, `rho_check.py`, `Wcoeffs.py` - focused diagnostics.