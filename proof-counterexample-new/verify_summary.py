"""Compact post-run audit of hunt outputs."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_lines(p):
    out = []
    if p.exists():
        for line in p.read_text(encoding="utf-8").splitlines():
            try:
                out.append(json.loads(line))
            except Exception:
                out.append({"parse_error": line[:80]})
    return out


def expected():
    spec = {}
    for j in (3, 4, 5):
        spec[j] = (801, 1200)
    for j in range(6, 61):
        spec[j] = (12 * j, 16 * j - 1)
    return spec


def main():
    spec = expected()
    rows = read_lines(HERE / "rows.jsonl")
    pairs = read_lines(HERE / "pairs.jsonl")
    susp = read_lines(HERE / "suspects.jsonl")
    tos = read_lines(HERE / "timeouts.jsonl")
    summary = json.loads((HERE / "summary.json").read_text(encoding="utf-8"))

    exp_rows = exp_pairs = 0
    by_j_rows, by_j_pairs = {}, {}
    problems = []
    seen_rows, seen_pairs = set(), set()
    for j, (lo, hi) in spec.items():
        exp_rows += hi - lo + 1
        exp_pairs += hi - lo + 1
    for r in rows:
        key = (r.get("j"), r.get("n"))
        if key in seen_rows:
            problems.append(f"duplicate row {key}")
        seen_rows.add(key)
        j, n = key
        if j not in spec or not (spec[j][0] <= n <= spec[j][1]):
            problems.append(f"out-of-spec row {key}")
        by_j_rows[j] = by_j_rows.get(j, 0) + 1
    for r in pairs:
        key = (r.get("j"), r.get("n"))
        if key in seen_pairs:
            problems.append(f"duplicate pair {key}")
        seen_pairs.add(key)
        j, n = key
        if j not in spec or not (spec[j][0] - 1 <= n <= spec[j][1] - 1):
            problems.append(f"out-of-spec pair {key}")
        by_j_pairs[j] = by_j_pairs.get(j, 0) + 1

    failures = [r for r in pairs if r.get("ok") is False and not r.get("skipped")]
    skipped = [r for r in pairs if r.get("skipped")]
    suspects = [r for r in susp if r.get("kind") == "row_suspect"]
    pair_fails = [r for r in susp if r.get("kind") == "pair_failure"]

    print(json.dumps({
        "expected_rows": exp_rows, "actual_rows": len(rows),
        "expected_pairs": exp_pairs, "actual_pairs": len(pairs),
        "missing_rows": sorted(set((j, n) for j, (lo, hi) in spec.items() for n in range(lo, hi + 1)) - seen_rows)[:20],
        "missing_pairs": sorted(set((j, n) for j, (lo, hi) in spec.items() for n in range(lo - 1, hi)) - seen_pairs)[:20],
        "failures": len(failures), "skipped_pairs": len(skipped),
        "suspects": len(suspects), "pair_failure_certificates": len(pair_fails),
        "timeouts": len(tos),
        "summary": summary,
        "problems": problems[:20],
    }, indent=2))


if __name__ == "__main__":
    main()
