#!/usr/bin/env python3
"""Independent exact definitions and certificates for ordered big-block rows.

Only the corrected recurrence is used:
  O_0 = 1,
  O_n = sum_{r=1}^{j-1} C(n,r) O_{n-r}
        + t sum_{r=j}^n C(n,r) O_{n-r}.

Coefficients are ascending in t.  The EGF implementation in egf_check.py is
independent of this recurrence.
"""
from __future__ import annotations

from contextlib import nullcontext
from dataclasses import dataclass
from functools import cmp_to_key
import math
from typing import Iterable, Iterator

from flint import acb, arb, ctx, fmpq, fmpq_poly, fmpz_poly


def ordered_rows(j: int, nmax: int) -> Iterator[tuple[int, list[int]]]:
    """Yield every corrected row O_{n,j} through nmax."""
    if j < 3:
        raise ValueError("this project deliberately skips j=2; require j>=3")
    rows: list[list[int]] = [[1]]
    yield 0, rows[0]
    for n in range(1, nmax + 1):
        deg = n // j
        row = [0] * (deg + 1)
        # Independently sum both parts of the stated binomial recurrence.
        for r in range(1, n + 1):
            c = math.comb(n, r)
            shift = int(r >= j)
            prev = rows[n - r]
            for m, a in enumerate(prev):
                row[m + shift] += c * a
        rows.append(row)
        yield n, row


def one_row(j: int, n: int, cache: dict[int, list[int]] | None = None) -> list[int]:
    """Convenience wrapper; a cache is needed to make one_row useful."""
    if cache is None:
        cache = {0: [1]}
    for k in range(len(cache), n + 1):
        row = [0] * (k // j + 1)
        for r in range(1, k + 1):
            c = math.comb(k, r)
            shift = int(r >= j)
            for m, a in enumerate(cache[k - r]):
                row[m + shift] += c * a
        cache[k] = row
    return cache[n]


def poly(coeffs: list[int]) -> fmpz_poly:
    return fmpz_poly(coeffs)


def squarefree_part(p: fmpz_poly) -> fmpz_poly:
    if p.degree() <= 0:
        return p
    q = p // p.gcd(p.derivative())
    # Positive scaling only: retain the Sturm sign convention.
    c = abs(q.content())
    return q // c if c else q


def primitive_positive_scale(p: fmpz_poly) -> fmpz_poly:
    c = abs(p.content())
    return p // c if c else p


def signed_pseudo_remainder(a: fmpz_poly, b: fmpz_poly) -> tuple[fmpz_poly, int]:
    """Return a positive multiple of -rem(a,b), and the pseudo-step count.

    Each elimination computes lc(b)*R - lc(R)*x^shift*b, hence after e
    eliminations R = lc(b)^e * rem(a,b), up to positive primitive scaling.
    The sign of lc(b)^e is therefore essential when forming the Sturm successor.
    """
    r = a
    lb = int(b.leading_coefficient())
    n = b.degree()
    e = 0
    while r.degree() >= n:
        lr = int(r.leading_coefficient())
        shift = r.degree() - n
        r = lb * r - lr * b.left_shift(shift)
        e += 1
        if r == 0:
            break
        r = primitive_positive_scale(r)
    if not r:
        return r, e
    if lb > 0 or e % 2 == 0:
        r = -r
    return primitive_positive_scale(r), e


@dataclass
class SturmChain:
    polynomials: list[fmpz_poly]
    degree: int
    squarefree_degree: int
    method: str = "pseudo"

    @classmethod
    def build(cls, p: fmpz_poly, method: str = "pseudo") -> "SturmChain":
        q = squarefree_part(p)
        if q.degree() <= 0:
            return cls([], 0, 0, method)
        seq = [q, q.derivative()]
        while seq[-1].degree() > 0:
            if method == "pseudo":
                rem, _ = signed_pseudo_remainder(seq[-2], seq[-1])
            elif method == "rational":
                rem = -(fmpq_poly(seq[-2]) % fmpq_poly(seq[-1]))
                if rem == 0:
                    break
                rem = primitive_positive_scale(rem.numer())
            else:
                raise ValueError(method)
            if rem == 0:
                break
            seq.append(rem)
        return cls(seq, p.degree(), q.degree(), method)

    @staticmethod
    def _at(q: fmpz_poly, x: int | None) -> int:
        # At x=None evaluate the sign at +infinity.
        if x is None:
            return 1 if q[q.degree()] > 0 else -1
        v = q(x)
        return (v > 0) - (v < 0)

    def variations(self, x: int | None = None) -> int:
        signs = [s for q in self.polynomials if (s := self._at(q, x))]
        return sum(a != b for a, b in zip(signs, signs[1:]))

    def minus_inf_variations(self) -> int:
        signs = [s for q in self.polynomials if (s := (1 if q[q.degree()] > 0 else -1) * (-1 if q.degree() % 2 else 1))]
        return sum(a != b for a, b in zip(signs, signs[1:]))

    def interval_counts(self, a: int, b: int) -> dict[str, int]:
        """Exact distinct squarefree roots in the indicated intervals.

        Endpoints are handled according to the usual Sturm convention; this
        routine is used at non-root endpoints in the validation pipeline.
        """
        va, vb = self.variations(a), self.variations(b)
        return {"closed": va - vb, "open_left": va - vb - int(self.polynomials[0](a) == 0),
                "open_right": va - vb - int(self.polynomials[0](b) == 0)}

    def summary(self) -> dict[str, int]:
        if not self.polynomials:
            return {"minus_inf": 0, "zero": 0, "plus_inf": 0,
                    "negative_distinct": 0, "zero_roots": 0,
                    "positive_distinct": 0, "real_distinct": 0,
                    "squarefree_degree": 0}
        vm = self.minus_inf_variations()
        v0 = self.variations(0)
        vp = self.variations(None)
        zero_roots = int(self.polynomials[0](0) == 0)
        return {
            "minus_inf": vm,
            "zero": v0,
            "plus_inf": vp,
            # V(-inf)-V(0) uses the right-closed Sturm convention, so remove
            # a root at zero to obtain the strictly negative count.
            "negative_distinct": vm - v0 - zero_roots,
            "zero_roots": zero_roots,
            "positive_distinct": v0 - vp,
            "real_distinct": vm - vp,
            "squarefree_degree": self.squarefree_degree,
        }


def ball_leq(a: arb, b: arb) -> bool:
    """Rigorous a <= b, false when the enclosures do not decide it."""
    return bool(a.upper() <= b.lower())


def compare_disjoint_real_balls(a: arb, b: arb) -> int:
    """Certified ordering; never falls back to floating-point midpoints."""
    if a.upper() < b.lower():
        return -1
    if b.upper() < a.lower():
        return 1
    raise ValueError("real-root balls are not rigorously disjoint")


def root_balls(coeffs: list[int], prec: int | None = None) -> tuple[fmpz_poly, list[arb], dict]:
    """Rigorous complex-root balls and projections for every root.

    An imaginary ball containing zero alone does NOT prove that its root is
    real. Pairwise-disjoint conjugate-invariant balls do (see below), as does
    an exact Sturm count. Once real-rootedness is certified, each real
    projection is a valid enclosure of that real root.
    """
    p = poly(coeffs)
    target_prec = int(prec) if prec is not None else ctx.prec
    if p.degree() <= 0:
        return p, [], {"prec_bits": target_prec, "imag_contains_zero_all": True,
                       "pairwise_disjoint": True, "all_real_by_disjoint_arb": True,
                       "root_order_certified": True,
                       "nonreal_enclosures": [], "distinct_roots": 0,
                       "roots_with_multiplicity": 0}
    work = ctx.workprec(target_prec) if prec is not None else nullcontext()
    with work:
        roots = p.complex_roots()
        bad = []
        projections = []
        balls = [z for z, _ in roots]
        for z, mult in roots:
            if not z.imag.contains(0):
                bad.append({"real": str(z.real), "imag": str(z.imag), "multiplicity": int(mult)})
            projections.append((z.real, int(mult)))
        pairwise_disjoint = all(not balls[i].overlaps(balls[k])
                                for i in range(len(balls)) for k in range(i))
        # If a nonreal root and its conjugate existed, their two disjoint balls
        # would both contain the same real part; but every imaginary interval is
        # assumed to contain zero. The balls would then overlap. Thus disjoint
        # enclosures + all imaginary intervals containing zero prove all-real.
        arb_all_real = (not bad) and pairwise_disjoint
        if arb_all_real:
            # Sort distinct-root entries by certified endpoint comparisons,
            # before expanding multiplicities. No float64 ordering is used.
            projections.sort(key=cmp_to_key(
                lambda u, v: compare_disjoint_real_balls(u[0], v[0])))
            vals = [a for a, mult in projections for _ in range(mult)]
        else:
            # No ordering is claimed for nonreal/undecided roots, and these
            # projections are never used to prove interlacing.
            vals = [a for a, mult in projections for _ in range(mult)]
        info = {"prec_bits": target_prec, "imag_contains_zero_all": not bad,
                "pairwise_disjoint": pairwise_disjoint,
                "all_real_by_disjoint_arb": arb_all_real,
                "root_order_certified": arb_all_real,
                "nonreal_enclosures": bad, "distinct_roots": len(roots),
                "roots_with_multiplicity": sum(m for _, m in roots)}
    return p, vals, info


def refine_nonreal_enclosures(coeffs: list[int], start_prec: int = 192,
                              max_prec: int = 8192) -> dict:
    """Refine until each root known nonreal by Sturm has imag excluding zero."""
    p = poly(coeffs)
    chain = SturmChain.build(p)
    counts = chain.summary()
    needed = p.degree() - counts["real_distinct"]
    at = int(start_prec)
    while True:
        _, _, info = root_balls(coeffs, at)
        found = sum(rec["multiplicity"] for rec in info["nonreal_enclosures"])
        if found >= needed:
            info.update({"ok": True, "expected_nonreal_multiplicity": needed,
                         "sturm": counts})
            return info
        if at >= max_prec:
            info.update({"ok": False, "expected_nonreal_multiplicity": needed,
                         "sturm": counts})
            return info
        at *= 2


def newton_check(coeffs: list[int]) -> dict:
    """Exact Newton inequalities, returning the minimum normalized slack."""
    d = len(coeffs) - 1
    nonpositive = [{"k": k, "coefficient": int(a)} for k, a in enumerate(coeffs) if a <= 0]
    if nonpositive:
        return {"ok": False, "nonpositive": nonpositive, "checks": 0}
    if d < 2:
        return {"ok": True, "nonpositive": [], "checks": 0}
    checks = []
    worst = None
    for k in range(1, d):
        lhs = int(coeffs[k]) ** 2 * k * (d - k)
        rhs = int(coeffs[k - 1]) * int(coeffs[k + 1]) * (k + 1) * (d - k + 1)
        ok = lhs >= rhs
        # Positive rational slack; no floating approximation is used.
        slack = fmpq(lhs - rhs, rhs)
        rec = {"k": k, "ok": ok, "lhs": lhs, "rhs": rhs}
        checks.append(rec)
        if not ok:
            return {"ok": False, "nonpositive": [], "checks": len(checks), "failure": rec}
        if worst is None or slack < worst[0]:
            worst = (slack, k)
    out = {"ok": True, "nonpositive": [], "checks": len(checks)}
    if worst is not None:
        out["minimum_slack"] = {"k": worst[1], "numerator": int(worst[0].numerator),
                                "denominator": int(worst[0].denominator)}
    return out


def check_real_rooted(coeffs: list[int], prec: int = 192, sturm: bool = True,
                      sturm_method: str = "pseudo", with_root_data: bool = False) -> dict:
    p, roots, info = root_balls(coeffs, prec)
    disc = int(p.discriminant()) if p.degree() > 0 else 1
    chain = SturmChain.build(p, method=sturm_method) if sturm else None
    arb_ok = bool(info["all_real_by_disjoint_arb"])
    if chain is None and not arb_ok:
        # Undecided/overlapping imaginary balls trigger an exact fallback;
        # contains(0) alone is never accepted as a real-root proof.
        chain = SturmChain.build(p, method=sturm_method)
    counts = chain.summary() if chain else None
    sturm_ok = (counts["real_distinct"] == chain.squarefree_degree) if chain else None
    # Two rigorous routes are available. Exact Sturm counting is preferred
    # when requested. Otherwise disjoint conjugate-invariant Arb root balls
    # prove all-real rootedness; merely containing zero in imag is not enough.
    # Discriminant >= 0 is an exact consistency cross-check in either route.
    certified = (sturm_ok is True) if chain is not None else arb_ok
    ok = certified and (disc >= 0)
    rec = {"degree": p.degree(), "ok": bool(ok), "arb": info,
           "discriminant_sign": (disc > 0) - (disc < 0),
           "discriminant_bits": abs(disc).bit_length(), "sturm": counts,
           "sturm_ok": sturm_ok, "arb_all_real_ok": arb_ok,
           "certification_method": ("exact_sturm" if sturm else
                                    "exact_sturm_fallback" if chain is not None else
                                    "disjoint_arb")}
    if not ok:
        rec["coefficients"] = [int(a) for a in coeffs]
        rec["nonreal_refinement"] = refine_nonreal_enclosures(coeffs, prec)
    rec["_discriminant"] = disc
    if with_root_data:
        rec["_root_data"] = (roots, info)
    return rec


def interval_interlacing(coeffs_f: list[int], coeffs_g: list[int], prec: int = 192,
                         max_prec: int = 4096, f_cert: dict | None = None,
                         g_cert: dict | None = None, f_root_data: tuple | None = None,
                         g_root_data: tuple | None = None) -> dict:
    """Rigorously compare consecutive ordered roots using Arb balls.

    Exact Sturm or disjoint-Arb certificates first establish that both rows
    are real-rooted. Strict root-ball inequalities then prove alternation.
    Equality/common-root cases
    are not guessed and are returned as inconclusive for separate treatment.
    """
    f, g = poly(coeffs_f), poly(coeffs_g)
    common = f.gcd(g)
    rec = {"degree_f": f.degree(), "degree_g": g.degree(),
           "common_degree": common.degree()}
    if g.degree() not in (f.degree(), f.degree() + 1):
        rec.update({"ok": False, "reason": "degree jump"})
        return rec
    fc = f_cert if f_cert is not None else check_real_rooted(coeffs_f, prec)
    gc = g_cert if g_cert is not None else check_real_rooted(coeffs_g, prec)
    rec["f_real_rooted"] = bool(fc["ok"])
    rec["g_real_rooted"] = bool(gc["ok"])
    if not (fc["ok"] and gc["ok"]):
        rec.update({"ok": False, "reason": "a row real-rootedness certificate failed"})
        return rec
    if common.degree() > 0:
        rec.update({"ok": None, "reason": "common roots; interval equality needs separate treatment"})
        return rec
    at = int(prec)
    while True:
        if at == int(prec) and f_root_data is not None:
            a, ia = f_root_data
        else:
            _, a, ia = root_balls(coeffs_f, at)
        if at == int(prec) and g_root_data is not None:
            b, ib = g_root_data
        else:
            _, b, ib = root_balls(coeffs_g, at)
        rec["arb_f"] = ia
        rec["arb_g"] = ib
        if len(a) != f.degree() or len(b) != g.degree():
            raise AssertionError("real-root count mismatch")
        if g.degree() == f.degree():
            # Equal degrees (1-based): f_1<=g_1<=f_2<=...<=f_d<=g_d.
            # The first loop includes the final f_d <= g_d comparison.
            comparisons = ([(a[i], b[i]) for i in range(f.degree())]
                           + [(b[i], a[i + 1]) for i in range(f.degree() - 1)])
        else:
            comparisons = [(b[i], a[i]) for i in range(f.degree())] + [(a[i], b[i + 1]) for i in range(f.degree())]
        rec["comparisons"] = len(comparisons)
        rec["certified"] = sum(ball_leq(x, y) for x, y in comparisons)
        if rec["certified"] == len(comparisons):
            rec.update({"ok": True, "prec_bits": at})
            return rec
        if at >= max_prec:
            rec.update({"ok": None, "reason": "insufficient root separation precision", "prec_bits": at})
            return rec
        at *= 2


def wronskian_crosscheck(coeffs_f: list[int], coeffs_g: list[int], prec: int = 192,
                         sturm: bool = True) -> dict:
    """Exact Wronskian plus Arb and exact Sturm real-root cross-check.

    This is deliberately NOT used as a stand-alone interlacing certificate.
    Certified ordered root-interval alternation is the primary check.
    """
    f, g = poly(coeffs_f), poly(coeffs_g)
    w = f * g.derivative() - f.derivative() * g
    rec = {"wronskian_degree": w.degree(), "common_degree": f.gcd(g).degree()}
    if w.degree() <= 0:
        rec["wronskian_ok"] = True
    else:
        wc = [int(a) for a in w]
        wc_rec = check_real_rooted(wc, prec=prec, sturm=sturm)
        rec.update({"wronskian_ok": wc_rec["ok"], "wronskian": wc_rec})
    rec["ok"] = bool(rec["wronskian_ok"])
    if not rec["ok"]:
        rec["coefficients_f"] = [int(a) for a in coeffs_f]
        rec["coefficients_g"] = [int(a) for a in coeffs_g]
    return rec
