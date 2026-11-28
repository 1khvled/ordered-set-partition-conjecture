# Mobius decoupling attempt (u = t/(1+t))

Date: 2026-11-28
Branch/work dir: `proof-mobius-decoupling`
Status: **obstruction located** (prompt identity false; correct coupled identity gives no clean Appell decoupling).

## 1. Setup

Fix $j\ge 2$. Let
$$
A_j(z)=\sum_{r=1}^{j-1}\frac{z^r}{r!},\qquad
C_j(z)=e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}=e^z-1-A_j(z),
$$
$$
F_j(z,t)=\frac1{1-A_j(z)-tC_j(z)}=\sum_{n\ge0}O_{n,j}(t)\frac{z^n}{n!}.
$$
The Mobius change is $u=t/(1+t)$, $t=u/(1-u)$, $1+t=1/(1-u)$.

## 2. Prompt identity is false; correct identity

The prompt claims $1-A_j-tC_j=(1+t)(1-A_j-u e^z)$. This is false. A direct simplification gives
$$
(1-A_j-tC_j) - (1+t)(1-A_j-u e^z)=2tA_j(z),
$$
since $(1+t)(1-A_j-u e^z)=(1+t)(1-A_j)-t e^z$ and
$1-A_j-tC_j=1-A_j-t(e^z-1-A_j)=(1+t)(1-A_j)-t e^z+2tA_j$.

Verified exactly by `verify.py` for $j=2..6$:
$j=2:2tz$, $j=3:tz(z+2)$, $j=4:tz(z^2+3z+6)/3$, etc., i.e. $2tA_j(z)$.

The correct identity is
$$
(1-u)(1-A_j(z)-tC_j(z)) = 1-A_j(z)(1-2u)-u e^z,
\qquad t=\frac{u}{1-u}. \tag{M}
$$
Equivalently
$$
F_j(z,t)=\frac{1-u}{1-A_j(z)(1-2u)-u e^z}. \tag{M'}
$$
Verified exactly in `verify.py` for $j=2..6$.

## 3. Correct coupled expansion

Let $B_j(u,z):=A_j(z)(1-2u)$. From (M'),
$$
F_j(z,t)=(1-u)\sum_{k\ge0}u^k \frac{e^{kz}}{(1-B_j(u,z))^{k+1}}. \tag{E}
$$
Using $(1-B_j)^{-(k+1)}=\sum_{p\ge0}\binom{k+p}{p}A_j(z)^p(1-2u)^p$, define
$$
 a_{m,p}^{(j)}:=m![z^m]A_j(z)^p\in\mathbb Z_{\ge0},
 \qquad a_{m,0}=\delta_{m,0},\ a_{m,p}=0\text{ unless }p\le m\le p(j-1),
$$
and
$$
 B_{n,k}^{(j)}(u):=n![z^n]\frac{e^{kz}}{(1-A_j(1-2u))^{k+1}}
 =\sum_{m=0}^{n}\binom{n}{m}k^{n-m}\sum_{p=0}^{m}\binom{k+p}{p}(1-2u)^p a_{m,p}^{(j)}. \tag{B}
$$
Then
$$
 O_{n,j}(t)=(1-u)\sum_{k\ge0}u^k B_{n,k}^{(j)}(u),\qquad t=u/(1-u), \tag{O}
$$
as formal power series in $u$ (the right side is rational/polynomial after cancellation).

## 4. Normalization and exact P formula

Let $q=\lfloor n/j\rfloor=\deg_t O_{n,j}$. The natural polynomial Mobius normalization is
$$
 P_{n,j}(u):=(1-u)^q O_{n,j}\!\left(\frac{u}{1-u}\right)\in\mathbb R[u],\qquad \deg P_{n,j}\le q.
$$
(The prompt wrote $(1+t)^q O_{n,j}(t)=(1-u)^{-q}O_{n,j}(u/(1-u))$, which is the reciprocal rational normalization; the polynomial one is $(1-u)^q O_{n,j}(u/(1-u))$.)

Combining with (O) gives the exact formula
$$
 \boxed{\;
 P_{n,j}(u)=(1-u)^{q+1}\sum_{k\ge0}u^k
 \sum_{m=0}^{n}\binom{n}{m}k^{n-m}\sum_{p=0}^{m}
 \binom{k+p}{p}(1-2u)^p\, m![z^m]A_j(z)^p\; } \tag{P}
$$
with the convention that the infinite $k$-sum collapses to a polynomial of degree $\le q$.

An equivalent coefficient form is obtained by expanding $(1-2u)^p=\sum_{i=0}^{p}\binom{p}{i}(-2)^i u^i$, $(1-u)^{q+1}=\sum_{s=0}^{q+1}\binom{q+1}{s}(-1)^s u^s$, and using Eulerian/falling-basis summation for $\sum_{k\ge0}k^r\binom{k+p}{p}u^k$. The resulting $u^r$ coefficient is an explicit finite integer combination of $a_{m,p}^{(j)}$, but it retains the coupled $(p,i)$ summation; no separation into a $k$-only Appell shift appears.

Verified by `verify.py` (`check_correct_decoupling`): (P) reproduces $P_{n,j}$ exactly for $j=3,4$, $n\le7$.

## 5. Root/interlacing under the Mobius map

For $t=u/(1-u)$, $u=t/(1+t)$:
- $t\in(-\infty,-1)\mapsto u\in(1,\infty)$;
- $t=-1$ is a pole ($u=\infty$);
- $t\in(-1,0]\mapsto u\in(-\infty,0]$.

Since $O_{n,j}$ has nonnegative coefficients and degree $q$, its real roots lie in $(-\infty,0]$. Hence $O_{n,j}$ real-rooted is equivalent to $P_{n,j}$ real-rooted with all real roots in $(-\infty,0]\cup(1,\infty)$ (no root in $(0,1)$ except possible endpoint behavior). Moreover strict interlacing of $O_n,O_{n+1}$ on $(-\infty,0]$ transfers to strict interlacing of $P_n,P_{n+1}$ on the corresponding real Mobius image, with the usual order reversal on intervals where the map is decreasing/increasing taken into account ($u=t/(1+t)$ is increasing on each component of $\mathbb R\setminus\{-1\}$).

This map is therefore equivalent, not simplifying: it preserves real-rootedness/interlacing as a statement.

## 6. Why the hoped decoupling fails

### 6.1 Algebraic obstruction

The prompt's clean decoupling requires the scalar-in-$z$ factorization
$$
1-A_j-tC_j = S(t)\,(1-A_j-u e^z).
$$
Comparing $e^z$-coefficients forces $S(t)=1+t$, but then the $A_j$-coefficient gives
$S(t)(1-A_j)=1-A_j+t(1+A_j)$, i.e. $2tA_j=0$, impossible for $j\ge2$. The exact residual is $2tA_j(z)$ (Section 2). Therefore there is no scalar Mobius factorization that isolates $e^{kz}/(1-A_j)^{k+1}$.

### 6.2 Coupling obstruction in the correct formula

The correct (M') does isolate $e^{kz}$, but over the $u$-dependent base
$1-B_j(u,z)=1-A_j(z)(1-2u)$:
$$
(1-A_j(1-2u))^{-(k+1)}=\sum_{p\ge0}\binom{k+p}{p}(1-2u)^pA_j(z)^p.
$$
Consequently $B_{n,k}^{(j)}(u)$ is not a polynomial in $k$ alone; it is a polynomial in $(k,u)$ with the coupling factor $(1-2u)^p$. The Appell-shift structure the prompt targets is exactly
$$
\sum_{k\ge0}u^k\sum_{m=0}^{n}\binom{n}{m}k^{n-m}d_m^{(k+1)},\qquad d_m^{(k+1)}=m![z^m](1-A_j)^{-(k+1)},
$$
which would follow only from the false identity. In the true identity, $d_m^{(k+1)}$ is replaced by $\sum_p \binom{k+p}{p}(1-2u)^p a_{m,p}^{(j)}$, re-entangling $u$ with the combinatorial $p$-sum.

### 6.3 No new leverage

After (P), proving real-rootedness/interlacing of $P_{n,j}$ is equivalent to the original conjecture (Section 5). The coefficient formula (P) does not remove the hard part: it moves the $t$-weight on the large blocks into a $(1-2u)^p$ weight on powers of $A_j$. Any positivity/interlacing argument in $u$ must still control the same block decomposition, now through binomial coefficients in $(k,p)$ and powers of $(1-2u)$. In particular, the Eulerian/Appell multiplier-sequence mechanism that would apply if $B_{n,k}$ depended only on $k$ is unavailable because of the $(1-2u)^p$ factor.

Thus the Mobius decoupling, as corrected, is an exact algebraic reformulation but not a structural simplification: it collapses to the same coupled block-counting inequality in a different basis.

## 7. Files

- `verify.py` - exact SymPy checks: prompt identity residual, correct identity, wrong decoupling mismatch, correct coupled formula reproduction.
- `NOTES.md` - this note.

## 8. Later correction: the decoupled variable is v = t/(1-t)

The parameter $u=t/(1+t)$ used above is the wrong Mobius map. Since
$C_j=e^z-1-A_j$,
\[
1-A_j-tC_j=(1+t)-(1-t)A_j-t e^z.
\]
Setting $v=t/(1-t)$, equivalently $t=v/(1+v)$, gives the exact decoupled form
\[
F_j(z,t)=\frac{1+v}{(1-A_j(z))-v(e^z-2)}
       =\frac{1+v}{1-A_j(z)-v e^z+2v}.
\]
Here $e^z-2$ is $j$-independent. The map $t\mapsto t/(1-t)$ is strictly
increasing on $(-\infty,0]$ and maps it onto $(-1,0]$, so root order and
consecutive interlacing are preserved. This corrects the earlier variable and
should be used instead of $u$ if the route is revisited. It still does not
constitute a proof.
