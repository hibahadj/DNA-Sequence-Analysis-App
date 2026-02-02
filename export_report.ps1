param(
    [string]$InFile = "REPORT.md",
    [string]$OutFile = "REPORT.docx"
)

$ErrorActionPreference = "Stop"

# Run from project root (folder containing this script)
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root

if (-not (Test-Path $InFile)) {
    Write-Error "Input file not found: $InFile"
}

$pandoc = Get-Command pandoc -ErrorAction SilentlyContinue
$pandocExe = $null
if ($pandoc) {
    $pandocExe = $pandoc.Source
} else {
    # Pandoc is often installed but not yet available in the current terminal PATH.
    $candidates = @(
        "C:\\Program Files\\Pandoc\\pandoc.exe",
        "C:\\Program Files (x86)\\Pandoc\\pandoc.exe",
        "C:\\Users\\$env:USERNAME\\AppData\\Local\\Pandoc\\pandoc.exe"
    )
    foreach ($p in $candidates) {
        if (Test-Path $p) {
            $pandocExe = $p
            break
        }
    }
}

if (-not $pandocExe) {
    Write-Host "Pandoc was not found on PATH or in common install locations." -ForegroundColor Yellow
    Write-Host "Install it with:" -ForegroundColor Yellow
    Write-Host "  winget install --id JohnMacFarlane.Pandoc -e" -ForegroundColor Yellow
    Write-Host "Then open a new terminal and re-run:" -ForegroundColor Yellow
    Write-Host "  .\\export_report.ps1" -ForegroundColor Yellow
    exit 2
}

& $pandocExe $InFile `
    --from=gfm `
    --to=docx `
    --standalone `
    --toc `
    -o $OutFile

Write-Host "Created $OutFile from $InFile" -ForegroundColor Green
