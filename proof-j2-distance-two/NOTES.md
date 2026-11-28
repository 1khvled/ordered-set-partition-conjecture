# `j=2` distance-two diagonal step: exact reduction and exact obstruction

Date: 2026-11-28  
Scope: only `ordered-big-blocks/proof-j2-distance-two` was created/edited.  
No installs, builds, external agents, models, or browser use. SymPy (already
installed) is used only for exact rational/Sturm verification.

## 0. Outcome

**No all-index proof of (I2) was obtained.** What is obtained is a complete,
rigorous determination of *why* the established toolkit cannot close it, and a
clean exact reduction of (I2) to a single new invariant:

\[
\boxed{
W(Q_{m,a},Q_{m-2,a+2})
=\frac{1}{a+1}\Bigl[\,W(Q_{m,a},\mathcal UQ_{m-2,a+1})
+a\,W(Q_{m,a},Q_{m-2,a+1})\Bigr].
}\tag{R}
\]

Identity (R) is a formal identity (residual zero for arbitrary polynomials,
§3). Hence (I2) is implied by the **distance-two coupled interlacing**

\[
\boxed{
W(Q_{m,a},Q_{m-2,a+1})>0
\quad\text{and}\quad
W(Q_{m,a},\mathcal UQ_{m-2,a+1})>0 .
}\tag{D2C}
\]

Rigorous result of this directory: the full conjunction of every fact proved
so far — (H), coefficient positivity, factorial-row real-rootedness, vertical
interlacing, the PF multipliers \(\Lambda_a\), the horizontal *coupled*
interlacing \(Q_{m-1,a}\preceq Q_{m,a}\), \(\mathcal UQ_{m-1,a}\preceq Q_{m,a}\),
the adjacent diagonal step (I1), and the vertical / \(\mathcal U\)-identities —
**does not imply (D2C), and therefore does not imply (I2)**. An exact rational
model satisfying all of these violates both premises of (D2C) and violates
(I2) (§5). This is an exact obstruction to the generic route, not a numerical
search.

The associated-Stirling recurrence is essential and cannot be replaced by any
of the above formal properties. The natural sign-count induction that proved
the adjacent coupled theorem fails at distance two for a specific, identified
reason (§6). Finite evidence through \(m=29\) (all \(a\)) finds no failure of
(D2C) or (I2) on the actual triangle; this is labelled finite evidence, never
a proof.

---

## 1. Setup and exact identities

Write \(A_{m,r}={m\brace r}_{\ge2}\) and

\[
Q_{m,a}(t)=\sum_{r\ge1}(a)_r A_{m,r}\,t^{r-1},\qquad
\mathcal U f:=2f+t f',\qquad W(f,g)=f'g-fg'.
\]

Three exact identities are used; all have zero symbolic residual (§7, V1–V3).

**(AS) Associated-Stirling recurrence.**
\[
A_{m+1,r}=rA_{m,r}+m A_{m-1,r-1}.
\tag{AS}
\]
Derived from \(\frac{d}{dz}\frac{C(z)^r}{r!}=\frac{C^r}{(r-1)!}+
\frac{zC^{r-1}}{(r-1)!}\) with \(C=e^z-1-z\), \(C'=C+z\). Verified exactly for
\(0\le m\le18\), \(0\le r\le8\): 171/171 residuals zero.

**(M) Polynomial recurrence.** From (AS), for \(m\ge2,a\ge1\):
\[
\boxed{Q_{m,a}=(tQ_{m-1,a})'+(m-1)a\,t\,Q_{m-2,a+1}.}
\tag{M}
\]
This is the recurrence that actually carries the associated-Stirling
structure; it is the recurrence used in `proof-u-operator/NOTES.md` §7 at
\(a=1\): \(G_m=(tG_{m-1})'+(m-1)t\,\mathcal UG_{m-2}\).

**(V) Vertical identity.**
\[
Q_{m,a+1}=Q_{m,a}+\frac1a(tQ_{m,a})'
=\frac1a\bigl(\mathcal UQ_{m,a}+(a-1)Q_{m,a}\bigr),
\tag{V}
\]
equivalently
\[
\boxed{\mathcal UQ_{m,a}=a\,Q_{m,a+1}-(a-1)Q_{m,a}.}
\tag{U}
\]
Both (V) and (U) have zero residual for \(4\le m\le15\), \(1\le a\le5\).

---

## 2. Reduction of (I2) to a coupled distance-two invariant

Put \(p=Q_{m,a}\), \(v=Q_{m-2,a+1}\), \(r=Q_{m-2,a+2}\). Applying (U) with
\(a\mapsto a+1\) to \(v\) gives
\[
\mathcal U v=(a+1)r-a v
\quad\Longleftrightarrow\quad
r=\frac{1}{a+1}\bigl(\mathcal U v+a v\bigr).
\]
Since \(W(p,\cdot)\) is linear in its second argument,
\[
W(p,r)=\frac{1}{a+1}\bigl[W(p,\mathcal Uv)+a\,W(p,v)\bigr],
\]
which is (R). Because \(a+1>0\) and \(a>0\), (R) has the **positive direction**
needed for §3 of the proof policy: if both brackets are positive, so is
\(W(p,r)\). Conversely, a failure of either bracket can still be masked by the
other, so (D2C) is sufficient but (R) shows it is the natural target, not an
equivalent one.

**Direction of (D2C).** In every nonautomatic case
\(\deg p=\deg v+1=\deg(\mathcal Uv)+1\). Hence the orientation corresponding
to \(W(p,v)>0\) and \(W(p,\mathcal Uv)>0\) is unambiguous: \(p\) surrounds \(v\)
and \(\mathcal Uv\). There is no sign ambiguity from equal-degree blocks at
distance two; equal degrees occur only between \(Q_{m,a}\) and \(Q_{m,a+1}\)
(same \(m\)), not along the diagonal step in question.

---

## 3. Symbolic residuals (policy item 5)

`verify.py` prints:

* (AS): 171/171 exact zero residuals.
* (V): 0 bad residuals for \(4\le m\le15\), \(1\le a\le5\).
* (U): 0 bad residuals over the same range.
* (R): residual exactly \(0\) for five random rational polynomial pairs and
  five random positive rational \(a\). Identity (R) is formal, not a coincidence
  of the associated-Stirling coefficients.

Finite enumeration is labelled finite evidence throughout; no finite check is
presented as a proof.

---

## 4. Direction and monotonicity audit (policy items 1–2)

**Direction samples (exact).** Representative outputs of `verify.py`:

| \((m,a)\) | \(W(p,v)\) | \(W(p,\mathcal Uv)\) | \(W(p,r)\) |
|---|---|---|---|
| \((7,1)\) | LC \(37800\), deg 2, 0 real roots | LC \(113400\), deg 2, 0 real roots | \(108(700t^2+35t+2)>0\) |
| \((8,1)\) | LC \(907200\), deg 4, 0 real roots | LC \(3628800\), deg 4, 0 real roots | \(18(126000t^4+84000t^3+38360t^2+880t+23)>0\) |
| \((9,2)\) | LC \(952560000\), deg 4, 0 real roots | LC \(3810240000\), deg 4, 0 real roots | \(16(119070000t^4+21168000t^3+2173290t^2+19866t+229)>0\) |
| \((10,3)\) | LC \(210039480000\), deg 6, 0 real roots | LC \(1050197400000\), deg 6, 0 real roots | \(450(933508800t^6+\cdots+43)>0\) |

All Wronskians have strictly positive leading coefficient and **zero** real
roots (exact Sturm count), hence are \(>0\) on all of \(\mathbb R\). The
direction is \(>0\), never \(<0\); (R) preserves it because the weights
\(1/(a+1)\) and \(a/(a+1)\) are positive.

**Discrete ratio monotonicity.** From (H) one gets the strict diagonal TP2
ratio chain \(R_0>R_1>R_2>\cdots\) with
\[
\frac{R_r}{R_{r-1}}
=\frac{a+r}{a+r-1}\frac{A_{m-1,r}A_{m,r-1}}{A_{m-1,r-1}A_{m,r}}
\le\frac{(a+r)(r-1)}{r(a+r-1)}<1\quad(r\ge2),
\]
and \(\frac{R_1}{R_0}\le\frac{a+2}{2(a+1)}<1\). Link-by-link direction is
printed in `proof-j2-diagonal-after-h/NOTES.md` §2; it is **not** used as a
substitute for (D2C).

**Caution — a reversed-direction trap.** For the *same* polynomial one has
\(W(Q_{m,a},\mathcal UQ_{m,a})=a\,W(Q_{m,a},Q_{m,a+1})<0\) (equal-degree
orientation; e.g. \(W(Q_{6,1},Q_{6,2})=-50-360t-4500t^2\)). This is the opposite
sign to (D2C) and to (I1). Any argument that silently transfers the sign of
\(W(p,\mathcal U p)\) to \(W(p,\mathcal U v)\) at distance two reverses the
direction and is invalid. This is recorded here as a retraction-style warning
in the sense of policy item 6: no claimed theorem in this directory rests on
that transfer.

---

## 5. Exact obstruction: feasible falsification (policy item 3)

Let \(\Sigma\) be the conjunction of **all** facts proved so far:

1. all active coefficients positive;
2. inequality (H) and the ratio chain of §4;
3. every \(Q_{m,a}\) real-rooted with negative roots; vertical interlacing
   \(Q_{m,a}\preceq Q_{m,a+1}\);
4. PF multipliers \(\Lambda_a\), \(Q_{m,a}=\Lambda_aG_m\), and \(\mathcal U\)-commutation;
5. horizontal coupled interlacing \(Q_{m-1,a}\preceq Q_{m,a}\) and
   \(\mathcal UQ_{m-1,a}\preceq Q_{m,a}\);
6. adjacent diagonal step (I1): \(W(Q_{m,a},Q_{m-1,a+1})>0\);
7. the polynomial identities (V) and (U) with zero residual.

The recurrence (AS)/(M) is deliberately **not** in \(\Sigma\); see §6.

**Theorem (exact obstruction).** \(\Sigma\not\models\) (D2C), and
\(\Sigma\not\models\) (I2).

*Witness.* The rational rows of `proof-j2-diagonal-after-h/NOTES.md` §5:
\[
A_{4,1}=1,\ A_{4,2}=\tfrac1{100};\quad
A_{5,1}=1,\ A_{5,2}=\tfrac3{10};\quad
A_{6,1}=1,\ A_{6,2}=\tfrac34,\ A_{6,3}=\tfrac1{20}.
\]
`verify.py` §V6 gives exactly
\[
D_1=1+\tfrac32t+\tfrac3{10}t^2,\qquad
D_2=2+\tfrac95t,\qquad
D_3=3+\tfrac3{25}t,
\]
\[
W(D_1,D_2)=\tfrac3{50}(9t^2+20t+20)>0,\qquad
W(D_2,D_3)=\tfrac{129}{25}>0,
\]
so (I1) holds. Yet
\[
W(D_1,D_3)=\tfrac3{250}(3t^2+150t+365),
\qquad \operatorname{disc}_t=\tfrac{8154}{3125}>0,
\]
so (I2) **fails** (two negative real zeros). More to the point of this
directory, the two premises of (D2C) also fail:
\[
W(p,v)=\tfrac3{500}(3t^2+200t+490),\quad
\operatorname{disc}_t=\tfrac{7677}{6250}>0,
\]
\[
W(p,\mathcal Uv)=\tfrac3{500}(9t^2+400t+970),\quad
\operatorname{disc}_t=\tfrac{28143}{6250}>0.
\]
Both discriminants are positive, so neither premise is \(>0\) on \(\mathbb R\).
Thus the obstruction hits (D2C) itself, not only its conclusion. The remaining
items of \(\Sigma\) for this witness — nonnegative (H) residuals
\((\tfrac7{25},0,\tfrac3{20},\tfrac3{100})\), real-rooted factorial rows, the
four coupled factorial-row Wronskians, unweighted horizontal pairs, the
\(Q_{6,a}\)-vs-\(Q_{5,a}\) pairs for \(a=1,2,3\), and zero (V)-residual — are
verified exactly in `proof-j2-diagonal-after-h/NOTES.md` §5; (V)/(U) hold for
*any* coefficient array, hence automatically here. \(\square\)

Consequences:

* (I1) + coupled horizontal interlacing + vertical identity does not compose
  to distance two; this is a second, sharper instance of the non-transitivity
  recorded in `proof-j2-resultants/NOTES.md` §4.3.
* The missing input is not another coefficient inequality of type (H) and not
  another root monotonicity: it is genuinely the distance-two coupling (D2C).

---

## 6. Why the associated-Stirling recurrence is essential, and the exact gap

The recurrence (M) is *not* part of \(\Sigma\); adding (AS) with its base cases
determines the triangle uniquely, so no finite feasible-model obstruction can
survive against it. Any successful proof must therefore use (M) essentially.
The only established mechanism that uses (M) is the sign-count induction of
`proof-u-operator/NOTES.md` §7. Here is exactly where it stops at distance two.

Let \(v=Q_{m-2,a+1}\), \(s=Q_{m-1,a}\), \(p=Q_{m,a}\), and let \(\lambda\) be a
root of \(\mathcal Uv\). By (M),
\[
p(\lambda)=(t s)'(\lambda).
\]
By (U), \((ts)'=\mathcal Us-s=a\bigl(Q_{m-1,a+1}-Q_{m-1,a}\bigr)\), so
\[
\boxed{\;p(\lambda)=a\bigl[Q_{m-1,a+1}(\lambda)-Q_{m-1,a}(\lambda)\bigr].\;}
\tag{GAP}
\]

**Correction from `proof-j2-qsign/NOTES.md` §1.** The displayed identity
\(p(\lambda)=(ts)'(\lambda)\) is false at roots of \(\mathcal Uv\): recurrence
(M) contains the extra term \((m-1)a\,t\,Q_{m-2,a+1}\), so actually
\(p(\lambda)=(ts)'(\lambda)+(m-1)a\,\lambda\,v(\lambda)\). The correct remaining
condition is the magnitude inequality (Sign) in `proof-j2-qsign/NOTES.md`, not
a sign of \(Q_{m-1,a}(\lambda)\) alone. The paragraphs below are retained as the
original, now superseded, diagnosis.

The coupled theorem supplies the sign of \(Q_{m-1,a+1}\) at \(\lambda\) (via
\(\mathcal UQ_{m-2,a+1}\preceq Q_{m-1,a+1}\)) and the sign of \(Q_{m-1,a}\) at
the roots of \(Q_{m-2,a+1}\) (via \(Q_{m-2,a+1}\preceq Q_{m-1,a}\)); but
\(\lambda\) is a root of \(\mathcal Uv\), not of \(v\). The sign of
\(Q_{m-1,a}(\lambda)\) is therefore **not controlled by any established
interlacing**, and (GAP) is a difference of two terms of individually known
but potentially cancelling sign. The §5 witness realises precisely this
cancellation: both premises of (D2C) fail there.

This is the exact gap. Closing it requires a new invariant that controls the
\(Q_{m-1,a}\)-sign at roots of \(\mathcal UQ_{m-2,a+1}\) — equivalently a direct
proof that \(\mathcal UQ_{m-2,a+1}\preceq Q_{m-1,a}\) *and* \(Q_{m-2,a+1}\preceq
Q_{m-1,a}\) hold with a common orientation strong enough to survive the
difference in (GAP). One natural candidate — the strict total positivity of
the matrix \((A_{m,r})\) — was checked for all \(1\times1,2\times2,3\times3\)
minors of the \(m=0..11,r=0..5\) block (`tp_probe.py`): zero negative minors.
That is finite evidence only; neither TP nor its implication for (D2C) is
proved here.

---

## 7. Boundary and degeneracy cases (policy item 4)

`verify.py` §Boundary prints:

* \(m=0,1\): \(Q_{m,a}=0\) — (I2) vacuous/omitted; the top term \(n!x^n\) is
  separate (`proof-j2-diagonal-after-h/NOTES.md` §1).
* \(m=2,3\): \(Q\) constant; Wronskian with a nonconstant partner is
  \(>0\) automatically, with a constant partner it is \(0\); omitted from strict
  (I2).
* \(m=4\): \(Q_{4,a}\) linear, \(Q_{2,a+2}\) constant, \(W=Q'_{4,a}\,Q_{2,a+2}>0\).
* \(a=1\): the \((a-1)\) term in (U) vanishes; (R) reduces to
  \(W(p,r)=W(p,\mathcal Uv)/2\), still positive-direction.
* Inactive \(A_{m,r}=0\) entries: still contribute the \(r=1\) term
  \((a)_1A_{m,1}=a\), so \(Q_{m,a}\neq0\) for \(m\ge2\); no active coefficient is
  used at a zero index.
* Equal degrees: none occur for the (I2) pair (degrees differ by 1), so the
  sign convention is unambiguous.

---

## 8. Finite evidence (explicitly not a proof)

`deepcheck.py 30` checks, with exact Sturm counts, all three of \(W(p,v)>0\),
\(W(p,\mathcal Uv)>0\), \(W(p,r)>0\) for \(6\le m\le29\), all relevant \(a\):
372/372 cases pass for each of the three statements, zero failures. This
extends the \(M\le18\) checks of `proof-j2-diagonal-after-h/NOTES.md` §8. It is
finite evidence only and must not be reported as a proof.

---

## 9. Conclusion

* (I2) is **reduced** to the distance-two coupled interlacing (D2C) by the exact
  formal identity (R), with positive direction.
* (D2C) — and hence (I2) — is **not implied** by the complete set of previously
  proved structural facts; an exact rational witness is given (§5).
* The associated-Stirling recurrence is essential; the precise point where the
  existing sign-count induction fails is the uncontrolled sign in (GAP) (§6).
* No all-index proof of (I2) was found. No counterexample on the actual
  triangle was found through \(m=29\).

A rigorous all-index proof would have to supply the new invariant identified
in (GAP) (or a TP/LGV proof of (D2C)); that is the exact remaining work.

---

## 10. Files and reproduction

* `NOTES.md` — this report.
* `explore.py` — exact construction of \(A_{m,r}\), \(Q_{m,a}\), Wronskians;
  first coupled-interlacing probe.
* `probe_u2.py` — probe of a *false* operator identity (kept to document the
  direction trap of §4); prints the nonzero residuals explicitly.
* `deepcheck.py` — deep finite Sturm check through \(m=29\).
* `tp_probe.py` — finite total-positivity probe of \((A_{m,r})\).
* `verify.py` — exact verification of (AS), (V), (U), (R), direction samples,
  the §5 witness, and boundary cases.

From this directory: `py -3.14 verify.py` (and `py -3.14 deepcheck.py 30`).
No packages installed, no builds, all files small.
