#!/usr/bin/env node
/**
 * Exact finite minor check for a_{n,m} = [t^m] O_{n,j}(t).
 *
 * This is evidence only, not a proof of total nonnegativity for all j,n.
 * It uses the first-block recurrence with binomial(n,r), never
 * binomial(n-1,r-1), and fraction-free Bareiss elimination over BigInt.
 */
"use strict";

function binom(n, r) {
  let v = 1n;
  for (let i = 0; i < r; i++) v = (v * BigInt(n - i)) / BigInt(i + 1);
  return v;
}

function determinant(input) {
  const a = input.map((row) => row.slice());
  const n = a.length;
  let previous = 1n;
  let sign = 1n;
  for (let k = 0; k < n - 1; k++) {
    let pivotRow = k;
    while (pivotRow < n && a[pivotRow][k] === 0n) pivotRow++;
    if (pivotRow === n) return 0n;
    if (pivotRow !== k) {
      [a[pivotRow], a[k]] = [a[k], a[pivotRow]];
      sign = -sign;
    }
    const pivot = a[k][k];
    for (let i = k + 1; i < n; i++) {
      const aik = a[i][k];
      for (let j = k + 1; j < n; j++) {
        a[i][j] = (a[i][j] * pivot - aik * a[k][j]) / previous;
      }
    }
    previous = pivot;
  }
  return sign * a[n - 1][n - 1];
}

function combinations(n, r) {
  const out = [];
  function visit(start, current) {
    if (current.length === r) {
      out.push(current);
      return;
    }
    for (let i = start; i <= n - (r - current.length); i++) {
      visit(i + 1, current.concat(i));
    }
  }
  visit(0, []);
  return out;
}

function rowsThrough(j, top) {
  const rows = [[1n]];
  for (let n = 1; n <= top; n++) {
    const row = Array(1 + Math.floor(n / j)).fill(0n);
    for (let r = 1; r <= n; r++) {
      const choose = binom(n, r);
      const shift = r >= j ? 1 : 0;
      for (let m = 0; m < rows[n - r].length; m++) {
        row[m + shift] += choose * rows[n - r][m];
      }
    }
    rows.push(row);
  }
  return rows;
}

/** Is there a perfect matching in the positive-support bipartite graph? */
function supportAdmissible(rowIndices, columnIndices, j) {
  const matched = Array(columnIndices.length).fill(-1);
  function augment(i, seen) {
    for (let q = 0; q < columnIndices.length; q++) {
      const m = columnIndices[q];
      if (seen[q] || (m !== 0 && rowIndices[i] < j * m)) continue;
      seen[q] = true;
      if (matched[q] < 0 || augment(matched[q], seen)) {
        matched[q] = i;
        return true;
      }
    }
    return false;
  }
  for (let i = 0; i < rowIndices.length; i++) {
    if (!augment(i, Array(columnIndices.length).fill(false))) return false;
  }
  return true;
}

function check(j, side = 10) {
  const rows = rowsThrough(j, side);
  const matrix = rows.map((row) =>
    row.concat(Array(side - row.length).fill(0n))
  );
  let tested = 0;
  let admissible = 0;
  let negative = null;
  let admissibleZero = null;

  outer:
  for (let size = 2; size <= side; size++) {
    const indices = combinations(side, size);
    for (const I of indices) {
      for (const J of indices) {
        tested++;
        const value = determinant(I.map((i) => J.map((k) => matrix[i][k])));
        if (value < 0n) {
          negative = { size, I, J, value: value.toString() };
          break outer;
        }
        if (supportAdmissible(I, J, j)) {
          admissible++;
          if (value === 0n && admissibleZero === null) {
            admissibleZero = { size, I, J };
          }
        }
      }
    }
  }
  return { j, side, tested, admissible, negative, admissibleZero };
}

for (const j of [2, 3, 4]) {
  console.log(JSON.stringify(check(j)));
}
