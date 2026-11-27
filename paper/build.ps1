$ErrorActionPreference = "Stop"
$python = "C:\Users\Abdelli\AppData\Local\Python\pythoncore-3.14-64\python.exe"
if (-not (Test-Path $python)) { $python = "python" }
& $python (Join-Path $PSScriptRoot "build_pdf.py")
