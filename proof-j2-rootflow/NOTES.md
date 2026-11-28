# j=2 root-flow repulsion: exact flow identities, direction audits, and the obstruction

Date: 2026-11-28
Scope: only `ordered-big-blocks/proof-j2-rootflow` was created/edited.
No installs, no builds. Every algebraic residual below is printed as `0` by `verify.py`.

## Outcome

**No all-`n` root-flow repulsion proof was obtained. The route is obstructed in
three precise ways.** The new rigorous content is:

1. An exact coefficient/residue identity for the normalized `x`-flow,
   `∂_x ψ_n = -u A_n` with `A_n` a polynomial with **strictly positive
   coefficients** (all `n`, all `x>0`), and the exact pairwise gap ODE in
   residue form.
2. A **direction error** in the stated Dyson condition: with `mu_k` the
   unnormalised `t`-roots, `d/dx(mu_{k+1}-mu_k)>0`, not `<0`; the condition
   holds only for the normalized roots `u_k=mu_k/x^2`.
3. An **exact Calogero–Moser obstruction**: the linear velocity representation
   `∂_x O_n=(a t+b)O_n'+c O_n` exists for `n<=5` and **does not exist for
   `n>=6`** (exact zero-dimensional solve, empty solution set).
4. A **monotonicity obstruction**: the residue function
   `f(u)=∂_xψ_n(u)/ψ_n'(u)` satisfies `f'<0` (Wronskian negative), not `f'>0`;
   repulsion comes from upward pole jumps at the critical points, and those
   jumps are not controlled by the sign of `f'`.
5. A **feasible-constraint falsification witness** satisfying every structural
   constraint that is actually proved here, for which the gap derivative is
   positive.

Finite `n` checks are labelled evidence, never a proof. The only all-`x`
case proved here is `n=4` (and `n<=3` is vacuous).

---

## 1. Exact setup and the flow identity

Write `C(z)=e^z-1-z` and

```
F(x,z,t) = 1/(1 - x z - t C(z)),      O_n(x,t) = n![z^n]F.
```

Introduce the **normalized** variables

```
w = x z,        t = x^2 u,
C_x(w) := x^2 C(w/x) = sum_{r>=2} w^r/(r! x^{r-2}),
G(x,w,u) := F(x, w/x, x^2 u) = 1/(1 - w - u C_x(w)),
psi_n(x,u) := x^{-n} O_n(x, x^2 u) = n![w^n]G.
```

`psi_n` has nonnegative coefficients in `u` and `1/x`, and
`psi_n -> Phi_n` (the Fibonacci-polynomial sequence of `proof-j2-base/NOTES.md`
§3) as `x -> infinity`. `psi_n(x,0)=n!>0`; hence any real roots are negative.

Exact derivative of the kernel:

```
d_x C_x(w) = (2x-w)e^{w/x} - 2x - w
           = - sum_{r>=3} (r-2) w^r/(r! x^{r-1}),          (E1)
d_x G = u (d_x C_x) G^2 ,      d_u G = C_x G^2 .           (E2)
```

Put `V_m(x,u) = m![w^m]G^2 = sum_{i=0}^m binom(m,i) psi_i psi_{m-i}`.
Extracting `n![w^n]` in (E2) and using (E1) gives the **exact flow identity**

```
    ∂_x psi_n(x,u) = -u A_n(x,u),
    A_n(x,u) := sum_{r=3}^n (r-2)/x^{r-1} * binom(n,r) * V_{n-r}(x,u).   (F)
```

Every coefficient of every `V_m`, hence of every `A_n`, is **nonnegative**
(and the `r=3` term makes `A_n` nonzero). `verify.py` prints residual `0` for
(F), for the scaling identity, and for `V_m = m![w^m]G^2` for `2<=n<=10`.

In the original `t` variable the same computation is
`∂_x O_n = n U_{n-1}`, `U_m = sum_i binom(m,i) O_i O_{m-i}`.

## 2. Exact pairwise gap ODE

Assume at a fixed `x>0` the polynomial `psi_n(x,.)` has `d=floor(n/2)` simple
roots `u_1<...<u_d<0`, and continue them as differentiable branches. Define

```
Q(x,u) := ∂_x psi_n(x,u) = -u A_n(x,u),
rho_k(x) := -u_k'(x) = Q(x,u_k)/psi_n'(x,u_k).               (G1)
```

The **residue identity** behind (G1) is

```
∂_x log psi_n(x,u) = (∂_x log lead) + sum_{k=1}^d rho_k/(u-u_k),   (G2)
```

so `rho_k` is exactly the residue of `∂_x log psi_n` at `u_k`. For adjacent
roots, with `g_k := u_{k+1}-u_k>0`,

```
   g_k' = u_{k+1}' - u_k' = rho_k - rho_{k+1}                      (G3)
        = [ Q(u_k) psi_n'(u_{k+1}) - Q(u_{k+1}) psi_n'(u_k) ]
          / [ psi_n'(u_k) psi_n'(u_{k+1}) ].                       (G4)
```

Because `psi_n'(u_k) psi_n'(u_{k+1})<0`, the required repulsion `g_k'<0`
is equivalent to **numerator of (G4) > 0**, i.e. to `rho_{k+1}>rho_k`.
Equivalently, substituting (F),

```
g_k' = [ u_k A_n(u_k) psi_n'(u_{k+1}) - u_{k+1} A_n(u_{k+1}) psi_n'(u_k) ]
       / [ psi_n'(u_k) psi_n'(u_{k+1}) ].                         (G5)
```

Equations (G3)–(G5) are the exact pairwise gap ODE requested. They are
identities; the problem is their sign.

## 3. Direction audit (policy item 1)

The task statement proposes `d/dx(mu_{k+1}-mu_k)<0` for unnormalised `t`-roots
`mu_k`. **This direction is false.** Since `mu_k(x)=x^2 u_k(x)`,

```
d/dx (mu_{k+1}-mu_k) = 2x g_k + x^2 g_k' .
```

At `x -> infinity`, `u_k(x) -> u_k(infty)<0` are distinct Fibonacci roots, so
`mu_{k+1}-mu_k ~ x^2 (u_{k+1}(infty)-u_k(infty))` with positive coefficient;
its derivative is therefore `>0` for all large `x`.

Explicit feasible sample (`n=4`, exact discriminant
`D=(36+8/x+1/x^2)^2-576`):

| `x` | `u`-gap | `d(u-gap)/dx` | `t`-gap | `d(t-gap)/dx` |
|---:|---:|---:|---:|---:|
| 1 | 6.344289 | **-1.970276** | 6.344289 | +10.718302 |
| 2 | 5.385326 | **-0.467126** | 21.541304 | +19.672801 |
| 10 | 4.651691 | **-0.018025** | 465.169055 | +91.231348 |

So the Dyson condition must be imposed on the **normalized** roots:
`d/dx(u_{k+1}-u_k)=g_k'<0`. With that correction, the condition is consistent
with the `x->infinity` data and is supported by finite evidence:
for `n=4,6,8,10` at `x=1,2,10` every sampled `g_k'` is negative
(`probe2.py`; e.g. `n=10,x=1`: `g'=(-15.09,-1.13,-0.23,-0.040)`).

## 4. Monotonicity audit (policy item 2)

Let `f(x,u):=Q(x,u)/psi_n'(x,u)`, so `rho_k=f(x,u_k)`. Then

```
f'(u) = N(u)/psi_n'(u)^2,     N := W(Q, psi_n') = Q' psi_n' - Q psi_n'' .   (M1)
```

At `x=1`, `verify.py` prints `-N` with **all positive coefficients** for
`3<=n<=10`; for `n=5,6` its discriminant is negative
(`-299555200`, `-3078999950400`), so `N<0` on all of `R`; e.g.

```
n=4: N=-450,      n=6: N=-20(59090 u^2 + 109618 u + 83405).
```

Thus `f` is **decreasing**, not increasing, between consecutive critical
points of `psi_n`. Let `c_k in (u_k,u_{k+1})` be the critical point. At `x=1`
the sample check gives `Q(c_k) psi_n''(c_k)>0` for every `k` (e.g. `n=4`,
`c_1=-15/4`, product `+450`), so `f` jumps from `-infinity` to `+infinity`
at `c_k`. Consequently

```
rho_{k+1}-rho_k = f(u_{k+1})-f(u_k)
```

is a competition between a **negative** integral of `N/psi_n'^2` on each side
of `c_k` and a **positive** pole jump. The pole term is
`~ 2 Q(c_k)/psi_n''(c_k) / eps` and cancels only in a Hadamard finite part;
the sign of `N` alone gives **no** bound on the finite remainder. So the
naive "increasing velocity function" sufficient condition is unavailable.

## 5. Exact Calogero–Moser obstruction at `n=6`

For `n<=5` there exist `a,b,c` (rational in `x`) with

```
∂_x O_n = (a t + b) O_n' + c O_n ,                            (CM)
```

so at a root `mu' = -(a mu+b)` and
`d(mu_{k+1}-mu_k)/dx = -a (mu_{k+1}-mu_k)`. For `n=4`, with
`den=(2x+1)(6x+1)(60x^2+8x+1)`,

```
a = -8(180x^3+108x^2+17x+1)/den,
b =  96x^3(4x+1)/den,
c =  16(180x^3+108x^2+17x+1)/den,
```

and `a<0`, so the `t`-gap grows: this agrees exactly with the direction audit
of §3 (it is not the `<0` condition stated in the task).

For `n>=6` the linear system for `(a,b,c)` is overdetermined (`deg_t O_n>=3`,
three unknowns, `>=4` equations) and `verify.py` prints an **empty solution
set** for `n=6,7,8,9,10` (residual-solve result `False`). Hence no representation
with velocity an affine function of the root exists from `n=6` onward. Identity
(F) is the correct replacement, but its velocity is a genuine
higher-order residue function.

## 6. Feasible-constraint falsification (policy item 3)

Everything proved above about a flow `Q=-uA` with `A` having positive
coefficients, negative simple roots, `N<0`, and upward pole jumps is shared by
the following exact witness:

```
psi(u) = (u+1)(u+2) = u^2+3u+2      (simple negative roots, positive coeffs),
Q(u)   = -u(50+31u)                  (A=50+31u has positive coefficients),
N(u)   = -2(31u^2+93u+75),  disc(N) = -2604 < 0  => N<0 on R,
critical point c=-3/2: Q(c) psi''(c) = 21/2 > 0  (upward jump).
```

Exact velocities:

```
rho_1 = Q(-2)/psi'(-2) = (-24)/(-1) = 24,
rho_2 = Q(-1)/psi'(-1) =   19 /  1  = 19,
g' = rho_1 - rho_2 = +5 > 0 .
```

So the gap **shrinks** as `x` decreases: the proposed repulsion fails. The
proved constraints are therefore **insufficient**; an extra family-specific
constraint (not identified here) is required. This is the exact obstruction.

## 7. Boundary and degeneracy cases (policy item 4)

* `n=0,1,2,3`: `d<=1`, there is no adjacent pair — vacuous.
* `n=4`: `d=2` and `A_4(x)=2(4x+1)/x^3>0`. With
  `s=36+8/x+1/x^2`, `D=s^2-576>0`,
  `rho_1-rho_2 = -A_4 s/(6 sqrt(D)) < 0` for **every** `x>0`.
  This is a complete all-`x` proof of repulsion for `n=4` (printed exactly by
  `verify.py` §5).
* `n=5`: `d=2`; the same degree-2 reduction applies but the all-`x` sign of
  `Q(r_1)+Q(r_2)` was not established here.
* `n>=6`: open; the CM form fails and the witness of §6 shows the available
  constraints do not close the gap.
* `x->infinity`: simple interlacing Fibonacci roots (imported from
  `proof-j2-base/NOTES.md` §3), so `g_k>0` initially. `x=1` is the endpoint of
  interest; no root collision is proved there.
* Coefficient zeros: the `r=2` term is absent from (F) because
  `d_x(x^{2-r})=0` at `r=2`; the `r=3` term is the leading positive term. No
  inactive index affects (F).

## 8. Symbolic residuals (policy item 5)

`verify.py` prints exact zero residuals for:
`x^{-n}O_n(x,x^2u)-psi_n`, `d_x psi_n + u*sum(...)`, and
`V_m - sum binom psi_i psi_{m-i}` for `2<=n<=10`; it prints empty solves for
the linear CM form at `n=6..10`; it prints the exact witness arithmetic.
All root computations in `probe2.py`/`probe4.py`/`probe6.py` are finite numerical
evidence only.

## 9. Retraction (policy item 6)

The stated sufficient condition `d/dx(mu_{k+1}-mu_k)<0` is retracted for the
original `t`-roots (direction error, §3). No theorem claiming all-`n`
repulsion is left in this record.

## 10. Files and reproduction

* `NOTES.md` — this report.
* `verify.py` — exact identity residuals, direction audit, CM obstruction,
  exact witness, `n=4` proof. Run `py -3.14 verify.py` here.
* `gen.py`, `probe1.py` ... `probe7.py` — exploration probes (finite evidence).
* `verify-out.txt` — saved output of `verify.py`.

No packages were installed and no builds were run.
