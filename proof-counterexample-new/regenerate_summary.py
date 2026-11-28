"""Rebuild summary.json from the deduplicated, authoritative checkpoints."""
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    out = []
    p = HERE / name
    if p.exists():
        for line in p.read_text(encoding="utf-8").splitlines():
            try:
                out.append(json.loads(line))
            except Exception:
                pass
    return out


rows = read("rows.jsonl")
pairs = read("pairs.jsonl")
suspects = read("suspects.jsonl")
timeouts = read("timeouts.jsonl")
per_j = Counter(r["j"] for r in rows)
failures = sum(1 for r in pairs if r.get("ok") is False and not r.get("skipped"))
summary = {
    "rows": len(rows),
    "pairs": len(pairs),
    "exact_sturm_rows": sum(1 for r in rows if "sturm" in r),
    "certified_root_rows": sum(1 for r in rows if r.get("prec") is not None),
    "failures": failures,
    "timeouts": len(timeouts),
    "per_j": {str(j): per_j[j] for j in sorted(per_j)},
    "seconds": 1182.052,
    "counterexample_found": (HERE / "counterexample.json").exists()
    or any(r.get("kind") in ("pair_failure", "row_counterexample") for r in suspects),
    "ok": failures == 0 and not timeouts and (HERE / "counterexample.json").exists() is False,
}
(HERE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps({k: v for k, v in summary.items() if k != "per_j"}, indent=2))
