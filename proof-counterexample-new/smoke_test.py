"""Offline smoke test for hunt.compute_job (uses a throwaway cache)."""
import json
import pickle
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import hunt
from scan import ordered_rows

tmp = HERE / "smoke_cache"
if tmp.exists():
    shutil.rmtree(tmp)
tmp.mkdir()
rows = [c for _, c in ordered_rows(6, 100)]
with (tmp / "j6.pkl").open("wb") as f:
    pickle.dump(rows, f)
# synthetic suspect rows: x^2+1 (disc<0, ULC fails)
with (tmp / "j7.pkl").open("wb") as f:
    pickle.dump([[1, 0, 1], [1, 1, 1], [1, 2, 1]], f)
hunt.CACHE = tmp

for j, n in [(6, 90), (6, 12), (7, 0), (7, 1)]:
    r = hunt.compute_job({"j": j, "n": n}, 64, 64)
    print(json.dumps({"j": j, "n": n,
                      "row": {k: v for k, v in r["row"].items() if k not in ("sha256",)},
                      "pair": r["pair"],
                      "sturm_used": r["sturm_used"],
                      "roots_used": r["roots_used"],
                      "counterexample": bool(r["counterexample"])}, indent=1), flush=True)
shutil.rmtree(tmp)
