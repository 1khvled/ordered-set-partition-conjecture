# Factorial-weighted 2-associated row: exact recurrence, Wronskian reduction, and the remaining lemma

Date: 2026-11-28
Scope: only `ordered-big-blocks/proof-factorial-row` was created/edited. No installs, builds,
other agents/models, or browser use.

## 0. Outcome

The corrected two-step differential recurrence for

$$E_r(x)=\sum_{k=1}^{d_r} k!\,A_{r,k}\,x^k,\qquad A_{r,k}=\Bigl\{{r\atop k}\Bigr\}_{\ge2},\quad d_r=\lfloor r/2\rfloor,$$

is verified exactly, and a **Wronskian criterion** reduces the desired real-rootedness of
$E_r$ to a sign inequality for a Wronskian $W_r$. I obtain an **exact decomposition**

$$W_r=B_1+(r-1)B_2$$

in which $B_1<0$ is proved unconditionally for all real $x$, and $B_2$ is reduced exactly to
one **operator lemma** (interlacing stability under $\mathcal U=2+x\,d/dx$). The operator lemma
is verified exactly for $r\le40$ and on 3000 random interlacing inputs, **but is not proved here**.
So this is **not yet an all-$r$ proof**; the obstruction is stated exactly in §6.

Everything in §§1–5 is rigorous; §7 records computational evidence only.

---

## 1. Exact recurrences (verified for $3\le r\le40$ in `explore.py`)

With $E_0=1,E_1=0$ and the three-term recurrence
$A_{r,k}=kA_{r-1,k}+(r-1)A_{r-2,k-1}$, multiplying by $k!x^k$ and summing gives

$$\boxed{\;E_r=xE_{r-1}'+(r-1)xE_{r-2}+(r-1)x^2E_{r-2}'\;}\tag{R}$$

equivalently $E_r=x\frac{d}{dx}\bigl[E_{r-1}+(r-1)xE_{r-2}\bigr]$.

Write $F_r:=E_r/x$ for $r\ge2$ (so $F_r(0)=1$, $\deg F_r=d_r-1$, all coefficients positive).
Then (R) is equivalent to each of

$$F_r=\frac{d}{dx}\Bigl[xF_{r-1}+(r-1)x^2F_{r-2}\Bigr],\tag{F}$$
$$U_r:=F_{r-1}+(r-1)xF_{r-2},\qquad F_r=(xU_r)'.\tag{U}$$

All three forms were checked exactly through $r=40$ (`explore.py`).

EGF cross-check: $\sum_{r\ge0}E_r(x)z^r/r!=\bigl[1-x(e^z-1-z)\bigr]^{-1}$, consistent with (R).

---

## 2. Wronskian criterion

_Definition._ $W(P,Q):=P'Q-PQ'$.

**Lemma 1.** If $W\bigl(F_r,\,xF_{r-1}\bigr)(x)<0$ for every real $x$, then $F_r$ and $xF_{r-1}$
have strictly interlacing real zeros; consequently $F_r$, hence $E_r$, has only real
non-positive zeros.

_Proof._ Put $R(x)=F_r(x)/(xF_{r-1}(x))$. Then $R'=W_r/(xF_{r-1})^2<0$ away from poles.
Near any pole $\alpha$, $R\sim c/(x-\alpha)$ and $R'<0$ forces $c>0$. Hence on each interval
between consecutive poles (the roots of $F_{r-1}$ and $0$), $R$ decreases from $+\infty$ to
$-\infty$ and has exactly one zero there; a degree count (using positive leading
coefficients) gives one extra zero of $F_r$ in $(-\infty,\alpha_1)$ exactly when
$\deg F_r=\deg(xF_{r-1})$, matching the observed parity pattern. Thus every zero of $F_r$ is
real and negative, and they alternate with the zeros of $xF_{r-1}$. $\blacksquare$

Once $E_r$ is real-rooted with non-positive roots, Newton's inequalities applied to
$c_k=k!A_{r,k}$ give (N') and hence (SLC) exactly as in
`proof-associated-stirling-lit/NOTES.md` §3 — the implication is already proved there and is
not repeated here.

---

## 3. Exact decomposition of $W_r$

Let $A:=F_{r-1}$, $C:=F_{r-2}$, so $F_r=A+xA'+(r-1)(2xC+x^2C')$. By (F),
$A=C+xC'+(r-2)(2xH+x^2H')$ with $H:=F_{r-3}$.
A direct expansion (verified symbolically for $5\le r\le40$ in `verify.py`) gives

$$\boxed{\;W_r=B_1+(r-1)B_2\;}\tag{1}$$
$$B_1:=x^2\bigl(AA''-A'^2\bigr)-A^2,\qquad
B_2:=x^2\bigl(3AC'-2A'C\bigr)+x^3\bigl(AC''-C'A'\bigr).$$

A second exact identity, valid for $r\ge7$ (verified for $7\le r\le40$), with
$\mathcal U f:=2f+xf'$, $D_C:=\mathcal U C$, $D_E:=\mathcal U E$ with $E:=F_{r-3}$, is

$$\boxed{\;B_2=x^2\,W(C,D_C)\;+\;(r-2)\,x^2\,W(D_C,\;xD_E)\;}\tag{2}$$

(the cases $r=5,6$ give $-B_2=12x^2$, $-B_2=22x^2$ respectively and are handled directly).

**Both identities are proved by elementary expansion of the definitions; no convergence or
induction hypothesis is used.**

---

## 4. The two manifestly negative pieces

### 4.1 $B_1<0$ everywhere (all $r$, no induction needed)

If $A$ has real roots $\alpha_i$ (all negative by the induction hypothesis), then
$A'/A=\sum_i\frac1{x-\alpha_i}$ and
$AA''-A'^2=-A^2\sum_i\frac1{(x-\alpha_i)^2}$. Hence

$$B_1=-A^2\left[1+x^2\sum_i\frac1{(x-\alpha_i)^2}\right]<0\qquad(x\in\mathbb R),$$

and at a root of $A$ one has $B_1=-x^2A'^2<0$ (note $A(0)=1$, so $x\ne0$ there).

### 4.2 $W(C,D_C)\le0$ everywhere

With $c=C'/C$ and $c'=-\sum_j\frac1{(x-\gamma_j)^2}$, and $D_C=2C+xC'=x^{-1}(x^2C)'$:

$$W(C,D_C)=C'(2C+xC')-C(2C'+C'+xC'')=-C^2\bigl(c-xc'\bigr).$$

Since $c-xc'=\sum_j\frac{-\gamma_j}{(x-\gamma_j)^2}>0$ for $C$ with negative roots
$\gamma_j$, we get $W(C,D_C)=-C^2\sum_j\frac{-\gamma_j}{(x-\gamma_j)^2}\le0$ for all real $x$,
with equality only at zeros of $C$. (Verified symbolically for $6\le r\le10$ in `verify.py`;
the derivation is general.)

---

## 5. Consequences already proved

Combining (1), (2), §4.1 and §4.2:

$$W_r=B_1+(r-1)\Bigl[x^2W(C,D_C)+(r-2)x^2W(D_C,xD_E)\Bigr],\qquad
B_1<0,\quad W(C,D_C)\le0.$$

Therefore

$$W_r<0\quad\text{for all real }x\quad\Longleftarrow\quad
\boxed{\;W\bigl(\mathcal U F_{r-2},\;x\,\mathcal U F_{r-3}\bigr)\le0\quad\text{for all real }x\;}\tag{L}$$

and by Lemma 1 this gives $E_r\in\mathrm{RZ}$ (non-positive roots), hence (N') and (SLC).

Note the induction structure: if the theorem holds for indices $<r$, then $A=F_{r-1}$ and
$C=F_{r-2}$ are real-rooted with negative roots, which is exactly what §4 needs. Thus the
**only** missing input for an all-$r$ proof is the single inequality (L), equivalently:

> **Operator lemma (missing).**
> Let $p=F_j$ and $q=F_{j-1}$ have strictly interlacing zeros with $xq$ (the orientation
> being: $p$ leftmost when $j$ is even, $xq$ leftmost and rightmost when $j$ is odd). Set
> $\mathcal U f=2f+xf'=(x^2f)'/x$. Then $\mathcal Up$ and $x\,\mathcal Uq$ have strictly
> interlacing negative zeros, equivalently
> $W(\mathcal Up,x\mathcal Uq)\le0$ for all real $x$.

The operator lemma is **false for arbitrary interlacing pairs of the wrong orientation**
(`op_test.py`: e.g. $p=x+43.6$, $q=x+16.3$), but appears true for exactly the two
orientations realised by this triangle (`op_test2.py`: 3000/3000 random cases).

---

## 6. Exact obstruction

* The real-rootedness of $E_r$ is **not proved** here.
* The reduction is complete: all that is missing is (L) / the operator lemma above.
* Equivalently, one may try to prove directly that
  $-W(\mathcal U F_{r-2},x\mathcal U F_{r-3})$ has no real zero; empirically it has all
  positive coefficients (`op_levels.py`, `w_positivity.py`) and exact Sturm counts confirm
  no real roots through $r=40$ (`verify.py`).
* Another route suggested by the data: the exact Wronskian ratio
  $$W(\mathcal U F_{r-2},\,x\mathcal U F_{r-3})
   =4\,W(F_{r-2},xF_{r-3})\,\frac{N_r(x)}{M_r(x)},\qquad
   -M_r=W(F_{r-2},xF_{r-3}),$$
  with $N_r,M_r$ having all positive coefficients and no real zeros (checked $r\le12$,
  `wq_ratio.py`, `w_factor.py`). Proving $N_r>0$ on $\mathbb R$ would also close the gap.
* A finite-descent structure exists: substituting $C=(xE)'+(r-3)(x^2H)'$ into
  $W(\mathcal UC,x\mathcal U E)$ produces a negative "self" term in $E$ plus a positive
  multiple of a Wronskian of the next operator level (`w_op_decomp.py`). Making this
  descent uniform in the level is the cleanest open path.

---

## 7. Computational evidence (not proof)

* `explore.py`: (R), (F), (U) verified exactly through $r=40$; $E_2,\dots,E_{10}$ printed.
* `verify.py`: identities (1), (2), §4.2 verified exactly through $r=40$; exact Sturm
  count confirms $-W_r$ has no real zero and $F_r$ alternates with $xF_{r-1}$ through
  $r=40$. This is the same finite range as `proof-associated-stirling-lit/NOTES.md` §4.
* `interlace_probe.py`, `structure_probe.py`, `pattern_probe.py`: strict alternation
  $F_r$ vs $xF_{r-1}$, and the fine pattern $u_i<\tau_i<\alpha_i<u_{i+1}<\dots$ where
  $u$ = zeros of $U_r$, $\tau$ = zeros of $F_r$, $\alpha$ = zeros of $F_{r-1}$.
* `op_test2.py`: the operator lemma holds for 3000/3000 random interlacing pairs in the
  two relevant orientations; `op_test.py` shows it fails in other orientations.
* `lemma_test.py`: naive "interlacing pencil + derivative" lemma is false (so Obreschkoff
  alone cannot close the step).

## 8. Files

`explore.py`, `verify.py`, `interlace_probe.py`, `structure_probe.py`, `wronskian_probe.py`,
`w_positivity.py`, `w_factor.py`, `w_identity.py`, `w_decomp_verify.py`, `w_op_decomp.py`,
`bracket_probe.py`, `b2_check.py`, `b2_identity.py`, `b2_coeffs.py`, `r_mono.py`,
`op_test.py`, `op_test2.py`, `op_levels.py`, `op_bracket.py`, `op_identity_search.py`,
`wq_ratio.py`, `wlog.py`, `wlog2.py`, `wlog3.py`, `lp_cert.py`, `n_cert.py`,
`n_recurrence.py`, `recurrence_search.py`, `identity_fit.py`, `d_recurrence.py`,
`pattern_probe.py`, `sos_search.py`, `wad.py`, `wop.py`, `lemma_test.py`.
