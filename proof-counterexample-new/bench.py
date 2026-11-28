import sys, time
from pathlib import Path
PARENT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PARENT))
from flint import ctx, fmpz_poly
from scan import ordered_rows, newton_failure

ctx.prec = 64
cases = [(3, 810), (4, 810), (5, 810), (6, 90), (30, 450), (60, 900)]
for j, target in cases:
    t0 = time.monotonic()
    gen = 0.0
    nf = 0.0
    disc = 0.0
    cnt = 0
    for n, coeffs in ordered_rows(j, target):
        t1 = time.monotonic()
        if n == target:
            cnt += 1
            t2 = time.monotonic()
            r = newton_failure(coeffs)
            t3 = time.monotonic()
            p = fmpz_poly(coeffs)
            d = p.degree()
            if 0 < d <= 64:
                dd = p.discriminant()
                s = (dd > 0) - (dd < 0)
            t4 = time.monotonic()
            gen += t2 - t1
            nf += t3 - t2
            disc += t4 - t3
            print(f"j={j} n={n} deg={d} newton={r and r.get('reason')} disc_sign={s if 0<d<=64 else 'skip'} gen_s={t2-t1:.3f} nf_s={t3-t2:.3f} disc_s={t4-t3:.3f}")
    print(f"j={j} TOTAL n=0..{target}: gen={gen:.3f}s nf={nf:.3f}s disc={disc:.3f}s wall={time.monotonic()-t0:.3f}s", flush=True)
