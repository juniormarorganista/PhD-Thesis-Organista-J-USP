# Windows Setup

This project contains many LaTeX files, images and long paths. On Windows,
this can cause failures during clone or extraction.

## 1. Enable long paths in Git

```powershell
git config --global core.longpaths true
```

## 2. Enable long paths in Windows

Open PowerShell as administrator:

```powershell
New-ItemProperty -Path "HKLM:\SYSTEM\CurrentControlSet\Control\FileSystem" `
  -Name "LongPathsEnabled" `
  -Value 1 `
  -PropertyType DWORD `
  -Force
```

Restart the computer.

## 3. Clone into a short path

Avoid OneDrive, Desktop and folders with long names.

Use:

```powershell
cd C:\
mkdir dev
cd C:\dev
git clone REPOSITORY_URL thesis
cd thesis
```

## 4. Install the Python layer

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements-dev.txt
python -m pip install -e .
python -m pytest -q
```

If virtual-environment activation is blocked:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
.\.venv\Scripts\Activate.ps1
```
