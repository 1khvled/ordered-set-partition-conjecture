# Middle-zero regime for Alexandersson's ordered big-block conjecture

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-middle` created/edited.  
Status: **no all-index middle-zero theorem proved**. A rigorous bulk saddle-point
framework is set up, the real-saddle large-powers asymptotic is identified, and
the precise obstruction to converting bulk asymptotics into uniform real-rootedness
is isolated.

## 1. Objects and normalizations

Fix `j>=2`, write `n=jm+s`, `0<=s<j`.  Let
```
A_j(z)=\sum_{r=1}^{j-1} z^r/r!,  C_j(z)=e^z-E_{j-1}(z),
G(z)=1/(1-A_j(z)),             H(z)=j! z^{-j} C_j(z) G(z).
```
Then `H(0)=1`, `H` has positive Taylor coefficients, and
```
c_{n,k} := [t^k]O_{n,j}(t)/n! = (j!)^{-k} [z^{n-jk}] H(z)^k G(z).
```
For the reversed normalized polynomial
```
Q_{m,s}(x)=\sum_{\ell=0}^{m} B_{m,\ell} (-x/D_m)^\ell,
B_{m,\ell}=c_{jm+s,m-\ell}/c_{jm+s,m},  D_m=D_j m^j, D_j=j!((j+2)/(j+1))^j,
```
the zeros `x` correspond to `t`-zeros by `t=-D_m/x`.  Left edge: fixed `ell`.
Right edge: fixed `m-ell`.  Middle regime: `ell=\alpha m`, `alpha\in(0,1)`.

## 2. Dominant-pole geometry of H

`H` is meromorphic. Its poles are zeros of `1-A_j(z)` (zeros of `C_j` do not
cancel generically; analytically H has the same pole set except possible removable
singularities, none observed).  For `j=2..6` the smallest-modulus zero of `1-A_j`
is the unique positive real zero `\alpha_j`, and all other zeros have strictly
larger modulus:

| j | alpha_j | next modulus |
|---|---:|---:|
| 2 | 1.000000 | - |
| 3 | 0.732051 | 2.732051 |
| 4 | 0.698885 | 2.930034 |
| 5 | 0.693903 | 3.254580 |
| 6 | 0.693233 | 3.591611 |

As `j->infty`, `alpha_j -> ln 2`.  The proof of uniqueness of the positive real
root is elementary: `A_j(0)=0`, `A_j` is strictly increasing on `[0,\infty)` with
`A_j(1)>1` for `j>=3` (`A_2(1)=1`).  The larger-modulus property is checked here
numerically for `j<=6` and is consistent with the limit `A_j -> e^z-1`; it is not
used as a theorem below except as motivation for the circle being inside the
dominant-pole disk.

## 3. Real-saddle large-powers setup

Let `ell=\alpha m`, `0<\alpha<1`, `k=m-ell=(1-\alpha)m`, `d=s+jell=s+j\alpha m`.
Then
```
I_{m,ell} := [z^d] H(z)^k G(z)
           = (1/2pi i) \oint e^{m F_\alpha(z)} G(z) z^{-s-1} dz,
F_\alpha(z) = (1-\alpha) log H(z) - j alpha log z.
```
Put
```
psi(z)=z H'(z)/H(z),      nu(z)=z psi'(z).
```
Since `H` has positive Taylor coefficients, for `z\in(0,\alpha_j)` the function
`psi` is the mean of a positive log-convex weight and is strictly increasing;
`nu(z)>0`.  Moreover `psi(0)=0` and `psi(z)->+infty` as `z->alpha_j^-` because
`H` has a simple pole at `alpha_j` (`A_j'(alpha_j)>0`).

**Lemma 3.1 (real saddle).** For every `alpha\in(0,1)` there is a unique
`z_alpha\in(0,alpha_j)` with
```
psi(z_alpha)=j alpha/(1-alpha)=:mu(alpha),   F_alpha'(z_alpha)=0.
```
*Proof.* `F_alpha'(z)=(1-\alpha)psi(z)/z - j alpha/z = (1-alpha)(psi(z)-mu(alpha))/z`.
Strict monotonicity and boundary values of `psi` give a unique root. `QED`.

Numerical values (bisection, `saddle_explore.py`) confirm the saddle moves from 0
to `alpha_j` as `alpha` increases.

## 4. Formal bulk saddle-point asymptotic

The standard large-powers/local-Gaussian formula (Hayman admissibility in the
positive-coefficient meromorphic setting; see e.g. Flajolet-Sedgewick, Analytic
Combinatorics, VIII.8, or Good 1957 for large powers) predicts, for fixed
`alpha\in(0,1)`,
```
[z^{d}]H(z)^k G(z) = H(z_alpha)^k G(z_alpha)
                     / (z_alpha^d sqrt(2pi k nu(z_alpha))) * (1+o(1)),
k=(1-alpha)m, d=j alpha m+s.
```
A numerical check (`saddle_asymptotic.py`, j=3, s=0) gives ratios actual/approx
close to 1 in the bulk and degrading at edges:

| m | ell | alpha | ratio |
|---:|---:|---:|---:|
| 20 | 5 | 0.25 | 0.991037 |
| 20 | 10 | 0.50 | 0.990749 |
| 20 | 15 | 0.75 | 0.983224 |
| 40 | 10 | 0.25 | 0.995513 |
| 40 | 20 | 0.50 | 0.995364 |
| 40 | 30 | 0.75 | 0.991569 |
| 40 | 39 | 0.975 | 0.922128 |

This supports the local Gaussian but does not prove it uniformly over all
`alpha\in(0,1)`.

## 5. Bulk rate function and failure of L_s normalization

The actual middle coefficient ratio is
```
B_{m,ell} = c_{n,m-ell}/c_{n,m}
          = (j!)^ell * ([z^{s+jell}]H^{m-ell}G / [z^s]H^m G).
```
The L_s normalization uses
```
B^L_{m,ell} = D_m^ell * s!/(s+jell)!,  D_m=D_j m^j.
```
For fixed `ell`, `B_{m,ell}/B^L_{m,ell}->1`.  For `ell=alpha m`, the saddle-point
prediction gives the bulk rate
```
R(alpha) := lim_{m->\infty} (1/m) log [z^{j alpha m}]H(z)^{(1-alpha)m}G(z)
          = (1-alpha) log H(z_alpha) - j alpha log z_alpha,
```
with `z_alpha\in(0,alpha_j)` real.  Since `H(z_alpha)>0` and `G(z_alpha)>0`, `R(alpha)` is
real.  The L_s normalization gives `(1/m)log B^L_{m,alpha m} = alpha j log m + O(1)`,
so the true bulk coefficients are exponentially larger in `m` than the L_s prediction,
matching `proof-final-tail/coeff_asympt.py`.  Therefore no fixed compact Hurwitz
argument based on `L_s` can control `ell=Theta(m)` zeros.

The real-valued rate means the formal bulk limit does **not** force complex zeros;
the obstruction is not a complex limit distribution but a uniformity/Szegő gap.

## 6. Why the bulk saddle-point does not give middle-zero real-rootedness

To turn the saddle-point into an all-index real-rootedness theorem one needs, in
addition to local Gaussian control, a global argument.  The standard routes and
their obstructions are:

1. **Uniform large-powers + zero-attracting real limit.** A uniform coefficient
   asymptotic `B_{m,alpha m} = exp(m R(alpha)+o(m))` with a smooth real rate `R`
   gives a limiting root distribution for the *scaled* polynomial, but not
   automatically real-rootedness of every finite `Q_{m,s}`.  The Szegő curve
   phenomenon is the canonical obstruction: partial sums of `e^z` have complex
   zeros despite an entire positive limit.  A real limit distribution is necessary
   but not sufficient.

2. **Circle-integral phase control.** Real-rootedness would follow if the Cauchy
   integral on `|z|=z_alpha` were dominated by the real saddle with a uniform
   error smaller than the saddle value, and if this dominance were compatible
   with a zero-free region in the relevant `x`-plane.  Proving this uniformly for
   all `alpha\in(0,1)`, including transitions to the edges, is exactly the missing
   global estimate.  The point `z_\alpha` approaches the pole `alpha_j` as
   `alpha->1` and approaches 0 as `alpha->0`; both ends are singular transition
   regimes.

3. **Riccati/linear-PDE/WKB route.** The linear PDE in
   `proof-final-inner-pde/NOTES.md` gives an exact recurrence but is not real
   stable and leaves a signed lagged derivative term.  Rewriting in the reversed
   `Q_{m,s}` variable produces a large-parameter ODE/PDE in `m` whose WKB turning
   points are precisely the bulk zeros.  Closing this would require an Airy/Bohr-
   Sommerfeld quantization with uniform error on the whole interval `alpha\in(0,1)`;
   the transition at `alpha->0` (left edge, L_s) and `alpha->1` (right edge,
   dominant-pole) are separate boundary layers.  No such uniform WKB theorem is
   proved here.

4. **Sign-regularity / Pólya frequency.** A PF or sign-regular kernel would give
   real-rootedness by standard closure theorems.  The kernel `H` is not PF in the
   required sense (see `literature.md` and prior proof notes); coefficient
   non-negativity alone is insufficient (Riccati counterexample in
   `proof/PROOF_NOTES.md` §7).  No sign-regularity structure closing the middle
   regime has been identified.

## 7. Precise obstruction

**Obstruction (middle-zero saddle-point gap).** Let `ell=alpha m` with
`alpha\in[\delta,1-\delta]`.  The local real saddle `z_alpha` exists and the
formal Gaussian approximation matches bulk coefficients numerically.  However, an
all-index real-rootedness theorem for `Q_{m,s}` requires one of the following
uniform controls, none proved here:

- (U1) A large-powers asymptotic for `[z^{s+jalpha m}]H^{(1-alpha)m}G` uniform in
  `alpha\in(0,1)`, with explicit error bounds up to the `alpha->0,1` boundary
  layers, and a global argument (e.g. Rouché on moving contours or a deformed
  saddle contour) excluding complex zeros of the full polynomial.
- (U2) A uniform WKB/Airy approximation for the reversed recurrence/ODE in `m`
  with matching boundary layers at the two edges and a Bohr-Sommerfeld condition
  forcing all roots real.
- (U3) A sign-regularity or variation-diminishing structure for the finite
  coefficient array that survives the `ell=Theta(m)` regime.

The edge theorems provide (U1)/(U2) only at fixed `ell` or fixed `m-ell`; the
Szegő phenomenon shows that local coefficient convergence is insufficient.  Thus
the middle-zero regime remains open after the saddle-point setup.

## 8. Reproduction

From `ordered-big-blocks/proof-final-middle`:

```powershell
py -3.14 .\saddle_explore.py
py -3.14 .\saddle_asymptotic.py
py -3.14 .\sturm_check.py
```

Files:
- `NOTES.md` - this note.
- `saddle_explore.py` - real saddle locations.
- `saddle_asymptotic.py` - local large-powers ratio check.
- `bulk_x_roots.py` - numpy roots (numerically unstable for large m; not a certificate).
- `bulk_x_roots_flint.py` - acb_poly roots (also unstable under this scaling; not a certificate).
- `sturm_check.py` - exact Sturm real-rootedness check for j=3, n<=120.
- `debug_psi.py`, `debug_bracket.py` - diagnostics.
