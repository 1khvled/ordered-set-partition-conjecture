"""Small exact obstruction probes for coefficient-matrix total positivity."""
from __future__ import annotations

from itertools import combinations
import json
from exact_algebra import ordered_rows


def det(matrix):
    n = len(matrix)
    if n == 0:
        return 1
    if n == 1:
        return matrix[0][0]
    # Fraction-free Bareiss elimination.
    a = [row[:] for row in matrix]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * pivot - a[i][k] * a[k][j]) // prev
            a[i][k] = 0
        prev = pivot
    return sign * a[-1][-1]


def main():
    found = []
    stats = {}
    for j in range(2, 7):
        rows = ordered_rows(j, 18)
        maxm = max(len(r) - 1 for r in rows)
        for k in range(2, min(5, maxm + 1) + 1):
            count = neg = zero = pos = 0
            for rr in combinations(range(len(rows)), k):
                for cc in combinations(range(maxm + 1), k):
                    matrix = [[rows[n][m] if m < len(rows[n]) else 0 for m in cc] for n in rr]
                    d = det(matrix)
                    count += 1
                    if d < 0:
                        neg += 1
                        if len(found) < 30:
                            found.append({"j": j, "k": k, "rows": rr, "cols": cc, "det": d})
                    elif d == 0:
                        zero += 1
                    else:
                        pos += 1
            stats[f"j={j},k={k}"] = {"count": count, "negative": neg, "zero": zero, "positive": pos}
    report = {"stats": stats, "first_negative_minors": found}
    from pathlib import Path
    Path(__file__).with_name("minor_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
