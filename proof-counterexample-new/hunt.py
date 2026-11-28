#!/usr/bin/env python3
"""Counterexample hunt driver.

All expensive root work runs in a killable worker subprocess so a single row
can be aborted by a HARD per-row timeout.  Cheap filters (Newton/ULC and
discriminant sign) decide which rows are suspects; Sturm and certified Arb
complex roots are used only for suspects.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

PARENT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PARENT))

from flint import ctx, fmpz_poly
from scan import newton_failure, sturm_variations, certified_complex_roots
from interlace import certainly_leq

OUTDIR = Path(__file__).resolve().parent
CACHE = OUTDIR / "cache"
WORKER = OUTDIR / "worker.py"
ROWS = OUTDIR / "rows.jsonl"
PAIRS = OUTDIR / "pairs.jsonl"
SUSP = OUTDIR / "suspects.jsonl"
TIMEOUTS = OUTDIR / "timeouts.jsonl"
SUMMARY = OUTDIR / "summary.json"
COUNTER = OUTDIR / "counterexample.json"
LOG = OUTDIR / "hunt.log"
WORKER_ERR = OUTDIR / "worker.err"


def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(line + "\n")


def digest(vals):
    return hashlib.sha256(json.dumps(vals, separators=(",", ":")).encode()).hexdigest()


def load_done(path):
    done = set()
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            try:
                r = json.loads(line)
                done.add((r.get("j"), r.get("n")))
            except Exception:
                pass
    return done


# --------------------------------------------------------------------------
# Computation (executed inside the worker process)
# --------------------------------------------------------------------------

_ROWS_CACHE = {}


def get_rows(j):
    """Load one cached row family; keep only the current one in memory."""
    if j not in _ROWS_CACHE:
        _ROWS_CACHE.clear()
        with (CACHE / f"j{j}.pkl").open("rb") as f:
            _ROWS_CACHE[j] = pickle.load(f)
    return _ROWS_CACHE[j]


def sturm_count(coeffs):
    vm = sturm_variations(coeffs, minus_inf=True)
    vp = sturm_variations(coeffs)
    return vm, vp, vm - vp


def real_root_enclosures(coeffs, prec0, max_prec=256):
    """Rigorous root balls; refine until non-reality/order is decided."""
    prec = prec0
    while prec <= max_prec:
        ctx.prec = prec
        roots = certified_complex_roots(coeffs, prec)
        nonreal = [(r, m) for r, m in roots if not r.imag.is_zero()]
        if nonreal:
            return prec, None, nonreal, False
        vals = [r.real for r, mult in roots for _ in range(mult)]
        # Sort by Arb midpoints, not float(): huge roots overflow float64.
        try:
            vals_sorted = sorted(vals, key=lambda a: a.mid())
        except Exception:
            vals_sorted = vals
        if all(certainly_leq(a, b) for a, b in zip(vals_sorted, vals_sorted[1:])):
            return prec, vals_sorted, [], True
        prec *= 2
    raise ArithmeticError(f"could not certify root order for degree {len(coeffs)-1}")


def interlace_from_roots(roots_f, roots_g):
    d, e = len(roots_f), len(roots_g)
    if e not in (d, d + 1):
        return False, {"ok": False, "reason": "degree_jump", "d": d, "e": e}
    ok = True
    if e == d:
        ok = all(certainly_leq(roots_f[i], roots_g[i]) for i in range(d))
        ok = ok and all(certainly_leq(roots_g[i], roots_f[i + 1]) for i in range(d - 1))
    else:
        ok = all(certainly_leq(roots_g[i], roots_f[i]) and certainly_leq(roots_f[i], roots_g[i + 1])
                 for i in range(d))
    rec = {"ok": bool(ok), "degree_f": d, "degree_g": e}
    if not ok:
        rec["bracket"] = {
            "f_roots": [str(r.mid()) for r in roots_f],
            "g_roots": [str(r.mid()) for r in roots_g],
        }
    return bool(ok), rec


def cheap_filter(coeffs, disc_limit):
    """Newton/ULC positivity test and exact discriminant sign.

    These are the ONLY computations applied to every row.
    """
    nf = newton_failure(coeffs)
    p = fmpz_poly(coeffs)
    degree = p.degree()
    disc_sign, disc_bits = None, None
    if 0 < degree <= disc_limit:
        try:
            disc = p.discriminant()
            disc_sign = (disc > 0) - (disc < 0)
            disc_bits = abs(int(disc)).bit_length()
        except Exception:
            disc_sign = "error"
    return {
        "nf": nf,
        "newton_ok": nf is None,
        "degree": degree,
        "disc_sign": disc_sign,
        "disc_bits": disc_bits,
        "suspect": (nf is not None) or (disc_sign == -1),
    }


def compute_job(job, prec, disc_limit):
    j, n, skip_pair = job["j"], job["n"], job.get("skip_pair", False)
    rows = get_rows(j)
    coeffs = rows[n]
    cf = cheap_filter(coeffs, disc_limit)
    degree, suspect = cf["degree"], cf["suspect"]

    # Expensive certificates: ONLY for suspects.
    roots = nonreal = prec_used = None
    root_error = None
    sturm_info = None
    if suspect:
        try:
            prec_used, roots, nonreal, _ = real_root_enclosures(coeffs, prec)
        except Exception as ex:
            root_error = f"{type(ex).__name__}: {ex}"
        try:
            vm, vp, cnt = sturm_count(coeffs)
            gcd_deg = p_gcd_degree(coeffs, degree)
            sturm_info = {
                "minus_inf": vm,
                "plus_inf": vp,
                "distinct_real": cnt,
                "gcd_with_derivative": gcd_deg,
                "complete": (cnt == degree) if gcd_deg == 0 else None,
            }
        except Exception as ex:
            sturm_info = {"error": f"{type(ex).__name__}: {ex}"}

    row = {
        "kind": "row",
        "j": j,
        "n": n,
        "degree": degree,
        "sha256": digest(coeffs),
        "newton_ok": cf["newton_ok"],
        "suspect": suspect,
        "prec": prec_used,
        "discriminant_sign": cf["disc_sign"],
        "discriminant_bits": cf["disc_bits"],
    }
    if not cf["newton_ok"]:
        row["newton_failure"] = cf["nf"]
    if nonreal:
        row["nonreal"] = [{"real": str(r.real), "imag": str(r.imag), "mult": m}
                          for r, m in nonreal]
    if root_error:
        row["root_error"] = root_error
    if sturm_info is not None:
        row["sturm"] = sturm_info

    suspect_rec = None
    if suspect:
        suspect_rec = {"kind": "row_suspect", "j": j, "n": n,
                       "coefficients": coeffs, "details": row}

    counterexample = None
    if nonreal or (sturm_info and sturm_info.get("gcd_with_derivative") == 0
                   and sturm_info.get("distinct_real") != degree):
        counterexample = {"kind": "row_counterexample", "j": j, "n": n,
                          "coefficients": coeffs, "details": row}

    pair = None
    if n > 0:
        if skip_pair:
            pair = {"ok": True, "j": j, "n": n - 1, "skipped": "prev_timeout"}
        else:
            prev_coeffs = rows[n - 1]
            prev_cf = cheap_filter(prev_coeffs, disc_limit)
            pair = check_pair(prev_coeffs, coeffs, prev_cf, cf, j, n,
                              prec, roots, nonreal, root_error)

    return {
        "j": j, "n": n, "row": row, "pair": pair,
        "suspect_rec": suspect_rec, "counterexample": counterexample,
        "sturm_used": sturm_info is not None,
        "roots_used": (roots is not None) or bool(nonreal) or (root_error is not None),
    }


def p_gcd_degree(coeffs, degree):
    p = fmpz_poly(coeffs)
    return p.gcd(p.derivative()).degree() if degree > 0 else 0


def check_pair(prev_coeffs, coeffs, prev_cf, cf, j, n, prec,
               roots, nonreal, root_error):
    do_interlace = prev_cf["suspect"] or cf["suspect"]
    if not do_interlace:
        return {"ok": True, "j": j, "n": n - 1, "skipped": True}
    prev_roots = prev_nonreal = prev_root_error = None
    if prev_cf["suspect"]:
        try:
            _, prev_roots, prev_nonreal, _ = real_root_enclosures(prev_coeffs, prec)
        except Exception as ex:
            prev_root_error = f"{type(ex).__name__}: {ex}"
    if roots is None and not nonreal and not root_error:
        try:
            _, roots, nonreal, _ = real_root_enclosures(coeffs, prec)
        except Exception as ex:
            root_error = f"{type(ex).__name__}: {ex}"
    if prev_nonreal or nonreal:
        ok, rec = False, {"ok": False, "reason": "nonreal"}
    elif prev_roots is None or roots is None:
        ok, rec = False, {"ok": False, "reason": "root_certification_unavailable",
                          "prev_error": prev_root_error, "curr_error": root_error}
    else:
        ok, rec = interlace_from_roots(prev_roots, roots)
    rec.update({"j": j, "n": n - 1, "prec": prec})
    return rec


# --------------------------------------------------------------------------
# Worker supervision
# --------------------------------------------------------------------------

def jobs_to_run(jmin, jmax, nmin=None, nmax=None):
    jobs = []
    if jmin <= 3 <= jmax:
        jobs.append((3, 801, 1200))
    if jmin <= 4 <= jmax:
        jobs.append((4, 801, 1200))
    if jmin <= 5 <= jmax:
        jobs.append((5, 801, 1200))
    for j in range(max(6, jmin), jmax + 1):
        jobs.append((j, 12 * j, 16 * j - 1))
    if nmin is not None or nmax is not None:
        clamped = []
        for j, lo, hi in jobs:
            lo2 = lo if nmin is None else max(lo, nmin)
            hi2 = hi if nmax is None else min(hi, nmax)
            if lo2 <= hi2:
                clamped.append((j, lo2, hi2))
        jobs = clamped
    return jobs


def reader(stream, q):
    try:
        for line in stream:
            q.put(line)
    except Exception:
        pass
    q.put(None)


def spawn_worker(args):
    err = WORKER_ERR.open("a", encoding="utf-8")
    proc = subprocess.Popen(
        [sys.executable, str(WORKER), "--prec", str(args.prec),
         "--disc-degree-limit", str(args.disc_degree_limit)],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err,
        text=True, bufsize=1,
    )
    q = queue.Queue()
    threading.Thread(target=reader, args=(proc.stdout, q), daemon=True).start()
    return proc, q


def kill_worker(proc):
    if proc is None:
        return
    try:
        proc.terminate()
        proc.wait(timeout=5)
    except Exception:
        try:
            proc.kill()
            proc.wait(timeout=5)
        except Exception:
            pass
    for stream in (proc.stdin, proc.stdout):
        try:
            stream.close()
        except Exception:
            pass


def ask(proc, q, job, cap):
    try:
        proc.stdin.write(json.dumps(job, separators=(",", ":")) + "\n")
        proc.stdin.flush()
    except Exception as ex:
        return None, f"send_failed: {ex}"
    try:
        line = q.get(timeout=cap)
    except queue.Empty:
        return None, "TIMEOUT"
    if not line:
        return None, "worker_eof"
    try:
        return json.loads(line), None
    except Exception as ex:
        return None, f"bad_result: {ex}; line={line[:200]!r}"


def write_line(path, rec):
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(rec, separators=(",", ":")) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prec", type=int, default=64)
    ap.add_argument("--disc-degree-limit", type=int, default=64)
    ap.add_argument("--row-time-cap", type=int, default=120)
    ap.add_argument("--jmin", type=int, default=3)
    ap.add_argument("--jmax", type=int, default=60)
    ap.add_argument("--nmin", type=int, default=None)
    ap.add_argument("--nmax", type=int, default=None)
    args = ap.parse_args()
    ctx.prec = args.prec

    jobs = jobs_to_run(args.jmin, args.jmax, args.nmin, args.nmax)
    done_rows = load_done(ROWS)
    done_pairs = load_done(PAIRS)
    start = time.monotonic()
    rows_count = pairs_count = sturm_rows = failures = timeouts = 0
    roots_rows = 0
    per_j = {}
    counterexample_found = False
    last_result = time.monotonic()
    current = {"j": None, "n": None}

    stop = threading.Event()

    def heartbeat():
        while not stop.wait(15):
            log(f"HEARTBEAT j={current['j']} n={current['n']} "
                f"since_last_result={time.monotonic()-last_result:.1f}s "
                f"rows={rows_count} pairs={pairs_count} elapsed={time.monotonic()-start:.1f}s")

    hb = threading.Thread(target=heartbeat, daemon=True)
    hb.start()

    proc, q = spawn_worker(args)
    log(f"START jmin={args.jmin} jmax={args.jmax} prec={args.prec} "
        f"row_cap={args.row_time_cap}s rows_done={len(done_rows)} pairs_done={len(done_pairs)}")

    try:
        for j, lo, hi in jobs:
            per_j[str(j)] = 0
            prev_timeout = False
            j_t0 = time.monotonic()
            for n in range(lo - 1, hi + 1):
                current.update(j=j, n=n)
                job = {"j": j, "n": n, "skip_pair": prev_timeout}
                result, why = ask(proc, q, job, args.row_time_cap)
                if result is None:
                    log(f"HARD TIMEOUT/FAIL j={j} n={n}: {why} "
                        f"elapsed={time.monotonic()-start:.1f}s")
                    kill_worker(proc)
                    with TIMEOUTS.open("a", encoding="utf-8") as f:
                        f.write(json.dumps({"j": j, "n": n, "reason": why,
                                            "seconds": round(time.monotonic() - start, 3)}) + "\n")
                    if n >= lo:
                        if (j, n) not in done_rows:
                            write_line(ROWS, {"kind": "row", "j": j, "n": n,
                                              "timeout": True, "reason": why})
                            done_rows.add((j, n))
                            rows_count += 1
                            per_j[str(j)] += 1
                        if (j, n - 1) not in done_pairs:
                            write_line(PAIRS, {"ok": True, "j": j, "n": n - 1,
                                               "skipped": "row_timeout"})
                            done_pairs.add((j, n - 1))
                            pairs_count += 1
                    timeouts += 1
                    prev_timeout = True
                    proc, q = spawn_worker(args)
                    last_result = time.monotonic()
                    continue

                prev_timeout = False
                last_result = time.monotonic()
                if result.get("error") or "row" not in result:
                    log(f"WORKER ERROR j={j} n={n}: {result.get('error', 'missing row')}")
                    if n >= lo:
                        write_line(SUSP, {"kind": "row_error", "j": j, "n": n,
                                          "error": result.get("error", "missing row")})
                    prev_timeout = True
                    continue
                if n >= lo:
                    if (j, n) not in done_rows:
                        write_line(ROWS, result["row"])
                        done_rows.add((j, n))
                        rows_count += 1
                        per_j[str(j)] += 1
                    if (j, n - 1) not in done_pairs and result["pair"] is not None:
                        rec = result["pair"]
                        write_line(PAIRS, rec)
                        done_pairs.add((j, n - 1))
                        pairs_count += 1
                        if not rec.get("ok") and not rec.get("skipped"):
                            failures += 1
                            write_line(SUSP, {"kind": "pair_failure", "j": j,
                                              "n": n - 1, "check": rec})
                            log(f"PAIR FAILURE j={j} n={n-1}: {rec.get('reason')}")
                    if result["suspect_rec"] is not None:
                        write_line(SUSP, result["suspect_rec"])
                        log(f"SUSPECT ROW j={j} n={n}")
                    if result["counterexample"] is not None and not counterexample_found:
                        counterexample_found = True
                        COUNTER.write_text(json.dumps(result["counterexample"], indent=2),
                                           encoding="utf-8")
                        log(f"COUNTEREXAMPLE j={j} n={n}")
                    if result["sturm_used"]:
                        sturm_rows += 1
                    if result["roots_used"]:
                        roots_rows += 1
                    log(f"j={j} n={n} rows={rows_count} pairs={pairs_count} "
                        f"sturm={sturm_rows} roots={roots_rows} failures={failures} "
                        f"timeouts={timeouts} elapsed={time.monotonic()-start:.1f}s")
            log(f"finished j={j} rows={per_j[str(j)]} total_rows={rows_count} "
                f"pairs={pairs_count} sturm={sturm_rows} roots={roots_rows} "
                f"failures={failures} timeouts={timeouts} "
                f"j_seconds={time.monotonic()-j_t0:.1f} elapsed={time.monotonic()-start:.1f}s")
    finally:
        stop.set()
        kill_worker(proc)

    summary = {
        "rows": rows_count,
        "pairs": pairs_count,
        "exact_sturm_rows": sturm_rows,
        "certified_root_rows": roots_rows,
        "failures": failures,
        "timeouts": timeouts,
        "per_j": per_j,
        "seconds": round(time.monotonic() - start, 3),
        "counterexample_found": counterexample_found,
        "ok": (failures == 0 and timeouts == 0 and not counterexample_found),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    log(f"SUMMARY {json.dumps(summary, separators=(',', ':'))}")


if __name__ == "__main__":
    main()
