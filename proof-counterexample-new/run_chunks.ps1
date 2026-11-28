$py = 'C:\Users\Abdelli\AppData\Local\Python\pythoncore-3.14-64\python.exe'
$script = 'proof-counterexample-new\hunt.py'
$log = 'proof-counterexample-new\chunks.log'
$chunk = 80

foreach ($j in 29..60) {
    $lo = 12 * $j
    $hi = 16 * $j - 1
    for ($start = $lo; $start -le $hi; $start += $chunk) {
        $end = [Math]::Min($start + $chunk - 1, $hi)
        & $py $script --jmin $j --jmax $j --nmin $start --nmax $end --row-time-cap 60 *>> $log
    }
}
& $py $script --jmin 3 --jmax 28 --row-time-cap 60 *>> $log
