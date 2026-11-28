"""Deduplicate rows/pairs checkpoints left by interleaved stale runs."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def clean(name):
    path = HERE / name
    seen = {}
    conflicts = []
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            rec = json.loads(line)
            key = (rec.get("j"), rec.get("n"))
        except Exception:
            continue
        if key in seen:
            if seen[key] != rec:
                conflicts.append(key)
            continue
        seen[key] = rec
    tmp = path.with_suffix(".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        for key in sorted(seen, key=lambda k: (k[0] if k[0] is not None else -1,
                                               k[1] if k[1] is not None else -1)):
            f.write(json.dumps(seen[key], separators=(",", ":")) + "\n")
    tmp.replace(path)
    return len(seen), conflicts[:10]


for name in ("rows.jsonl", "pairs.jsonl"):
    kept, conflicts = clean(name)
    print(f"{name}: kept={kept} conflicts={conflicts}")
