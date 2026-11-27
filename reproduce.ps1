param(
    [int]$SmallNMax = 600,
    [int]$JMax = 60
)

$ErrorActionPreference = "Stop"
$python = "C:\Users\Abdelli\AppData\Local\Python\pythoncore-3.14-64\python.exe"
if (-not (Test-Path $python)) { $python = "python" }

& $python -m pip install -r requirements.txt
& $python -m unittest discover -s tests -v
& $python scan.py --jmin 3 --jmax 20 --nmax 150 --out scan-results.jsonl
& $python interlace.py --jmin 3 --jmax 20 --nmax 150 --wronskian --out interlace-results.jsonl
& $python frontier.py --small-nmin 141 --small-nmax $SmallNMax --jmin 6 --jmax $JMax --out frontier-results.jsonl
