import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from scan import ordered_rows

# Generate rows for j=2..5
for j in [2,3,4]:
    rows = [c for _, c in ordered_rows(j, 40)]
    print('=== j=',j,'n rows:', len(rows))
    # print a few
    for n in range(0, 12):
        print(n, rows[n])
