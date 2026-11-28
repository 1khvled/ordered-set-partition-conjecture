# Proof-agent acceptance policy

Applies to every agent producing a claimed proof in this repository.

A result is **not accepted** merely because identities simplify to zero. Before
writing “proved”, the agent must complete all of the following and show the
evidence:

1. **Direction audit.** State every inequality with its exact direction and
   evaluate it at an explicit feasible sample point. A chain of inequalities
   must be checked one link at a time.

2. **Monotonicity audit.** If a function is claimed increasing or decreasing,
   print the derivative (or finite difference) and its sign on the claimed
   interval. Then check which endpoint gives the upper or lower bound. The
   C1 retraction was caused by reversing this direction.

3. **Feasible-constraint falsification.** Try to violate the claimed bound using
   only the constraints that were actually proved. If the constraints alone
   permit a counterexample, the argument is incomplete. Report the witness.

4. **Boundary and degeneracy cases.** Enumerate zero coefficients, equal
   degrees, top rows, and inactive indices. Show each is covered or excluded.

5. **Symbolic residuals.** Print zero residuals for every algebraic identity
   used. Finite enumeration is labelled finite evidence, never a proof.

6. **Retraction rule.** A sign error, reversed monotonicity, or unjustified
   optimization forces an immediate retraction in `NOTES.md` with the exact
   gap recorded. Do not leave an invalid “Theorem” in the record.

Violations are treated as failed work, not partial success.
