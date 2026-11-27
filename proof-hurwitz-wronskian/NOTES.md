# Hurwitz stability of the consecutive Wronskian W_n

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-hurwitz-wronskian`.  
Outcome: **no unconditional proof found; exact obstruction identified.**

## 1. Target

Fix `j>=2`. Let `O_n=O_{n,j}` be Alexandersson's ordered big-block polynomials.
Define the consecutive Wronskian

\[
W_n(t)=W(O_n,O_{n+1})(t)
=O_n(t)O_{n+1}'(t)-O_n'(t)O_{n+1}(t).
\]

Target (HW): for every `j>=2, n>=0`, every zero of `W_n` lies in the open
left half-plane `Re z<0`.

Finite exact checks (see `proof-automated-consecutive/NOTES.md`) show that for
tested ranges `W_n` has positive coefficients, no real zeros, and all complex
zeros have negative real parts.  These checks are evidence only.

## 2. General lemma: interlacing implies Hurwitz

There is a classical result of de Bruijn for Wronskians of real-rooted
polynomials.

**Lemma (de Bruijn).** Let `P,Q` be real-rooted polynomials with positive
leading coefficients whose roots interlace in the standard sense
(`deg Q in {deg P,deg P+1}`). Then `W(P,Q)=PQ'-P'Q` is a Hurwitz
polynomial: all of its zeros have strictly negative real parts.

*Proof sketch.* Write `R=Q/P`. Then `W(P,Q)=-P^2 R'`. The interlacing
hypothesis makes `R` a real rational function whose zeros and poles are real,
simple, and interlacing on the negative real axis. Such an `R` is a
Nevanlinna/Pick function with respect to the left half-plane; equivalently,
its critical points (the zeros of `R'`) are forced into `Re z<0`.  The
Hermite-Biehler formulation gives the same conclusion: the even and odd parts
of `W(P,Q)` have real negative roots which interlace.  QED.

Consequence for the family: if `O_n` and `O_{n+1}` are real-rooted and
consecutively interlacing, then `W_n` is automatically Hurwitz.

## 3. Reduction to the known open gaps

The all-parameter real-rootedness/interlacing of `O_n` is not proved.  The
available reductions are:

* **WRS (lag-`j` real-rootedness).** For `q=O_{n-j}`, `p=O_{n-1}`,
  real-rootedness of `O_n` is reduced by `proof-weighted/NOTES.md` to
  \[
  -\frac{R_n(p)(\mu)+S_n(\mu)}{q'(\mu)}
  >\binom{n-1}{j-1}\mu(\mu-1)
  \]
  at every root `\mu` of `q`.  This is the weighted root-separation inequality
  (WRS).

* **CONSEC-SIGN (consecutive interlacing).** By `proof-final-consecutive/NOTES.md`,
  if `O_n` is real-rooted with simple roots `\mu_i`, then
  \[
  O_n\interl O_{n+1}
  \iff
  \frac{O_{n+1}(\mu_i)}{O_n'(\mu_i)}<0
  \iff
  W_n(\mu_i)>0.
  \]
  This is the consecutive domination inequality (CONSEC-SIGN).

Thus (HW) for all `j,n` is implied by (WRS)+(CONSEC-SIGN), both of which
remain open.  The present target does not introduce a new independent gap; it
is the Hurwitz face of the same two missing quantitative inequalities.

## 4. Routh-Hurwitz route

For a polynomial with positive coefficients, Hurwitz stability is equivalent
to positivity of all leading principal minors of its Hurwitz matrix.  For
`W_n`, the coefficients are positive in all tested cases, and the Routh
array is positive in all tested cases.

However, the Routh minors are rational expressions in the coefficients of
`O_n` and `O_{n+1}`.  Their positivity is not a coefficient-ratio statement
and does not follow from the common compatible window: the same obstruction
that defeats WRS (the lagged `tQ` term, whose derivative part does not vanish
at the lagged row's zeros) reappears in the minors.  No recurrence with
polynomial-in-`n` coefficients for the Routh array was found, and no
stability-preserving recurrence for `W_n` alone was found.

## 5. Hermite-Biehler route

Write `W_n(t)=A_n(t^2)+tB_n(t^2)`.  By the Hermite-Biehler theorem,
`W_n` is Hurwitz iff `A_n` and `B_n` have only real negative roots and
their roots interlace.  Numerical checks show this for small cases.

Proving it in general again requires controlling the critical points of
`R=O_{n+1}/O_n`.  The even/odd root interlacing of `A_n,B_n` is not
accessible from the coefficient recurrence alone; it is governed by the same
CONSEC-SIGN domination.

## 6. Bivariate PDE / transfer function route

The exponential generating function is
\[
F_j(z,t)=\frac1{1-B_j(z,t)},\qquad
B_j(z,t)=\sum_{r=1}^{j-1}\frac{z^r}{r!}
+t\left(e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}\right).
\]
It satisfies the Riccati identity
\[
F_z=(2+(t-1)w_j)F^2-F,\qquad w_j=z^{j-1}/(j-1)!,
\]
and `F_t=C_jF^2`.  The coefficient `O_n` is `n![z^n]F`.

The Riccati equation integrates explicitly:
\[
F_j(z,t)=\frac{e^{-z}}{(2-t)+?}
\]
More precisely,
\[
\frac1{F_j(z,t)}=(t-1)E_{j-1}(z)-t e^z+2,
\qquad E_{j-1}(z)=\sum_{r=0}^{j-1}\frac{z^r}{r!}.
\]
This is a useful closed form, but it does not remove the obstruction: the
denominator is still linear in `t` with coefficients mixing `e^z` and the
truncated exponential `E_{j-1}`.  The factor `(t-1)` remains attached to the
truncation, and no coefficient extraction from this form gives a
stability-preserving recurrence for `W_n`.  Direct sampling also shows `F` is
not real stable in the ordinary sense, so the Borcea-Branden coefficient
theorem cannot be used without a nontrivial transformation.

## 7. Exact obstruction

The following exact obstruction is taken from `proof-weighted/NOTES.md`,
Section 6, and is the reason no compatibility-only proof can work.

Take `j=2`, `n=5`, so `c=4`, and
\[
p(t)=(t+2)(t+4)=t^2+6t+8,\qquad q(t)=t+3.
\]
The roots of `p` are `-4,-2`; the root `-3` of `q` lies between them, so
`(q,p)` is an ordinary compatible pair.  Let
\[
A_5=R_5(p)=5p+2tp'=9t^2+42t+40.
\]
Its roots are `-10/3,-4/3`, and `-10/3<-3<-4/3`; hence `q` is also
compatible with `A_5`.  Thus both ordinary compatibility and the proposed
transformed-window conclusion hold.

Nevertheless, at `\mu=-3`,
\[
p(\mu)=-1,\quad p'(\mu)=0,\quad q'(\mu)=1,\quad \alpha=1,\quad f=5,
\]
while
\[
c\mu(\mu-1)=4(-3)(-4)=48.
\]
So WRS fails: `5>48` is false.  The last term is nonzero:
\[
tQ(q)(-3)=(-3)(-4)(1)=12.
\]
The final combination is
\[
A_5+4tQ(q)=17t^2+50t+40,
\]
whose discriminant is `50^2-4\cdot17\cdot40=-220<0`.  Hence it has two
nonreal roots.

Conclusion of the obstruction: no theorem that uses only the common compatible
window, the structural coefficient `c`, and compatibility of `q` with `A_n`
can prove the required root step.  A quantitative invariant such as WRS or
CONSEC-SIGN is indispensable.  The same obstruction blocks the Routh-Hurwitz
and Hermite-Biehler routes for `W_n`, because those routes also need the
lagged derivative term to be controlled quantitatively.

## 8. Status

- (HW) is implied by (WRS)+(CONSEC-SIGN), both open.
- No unconditional proof of (HW) was obtained.
- No counterexample to (HW) was found in the tested ranges.
- Finite checks remain evidence only.
PY