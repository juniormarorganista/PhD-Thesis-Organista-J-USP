# Configuração no Windows

Este projeto contém muitos arquivos LaTeX, imagens e caminhos longos. No Windows, isso pode causar falhas no clone ou na extração.

## 1. Ativar caminhos longos no Git

```powershell
git config --global core.longpaths true
```

## 2. Ativar caminhos longos no Windows

Abra o PowerShell como administrador:

```powershell
New-ItemProperty -Path "HKLM:\SYSTEM\CurrentControlSet\Control\FileSystem" `
  -Name "LongPathsEnabled" `
  -Value 1 `
  -PropertyType DWORD `
  -Force
```

Reinicie o computador.

## 3. Clonar em caminho curto

Evite OneDrive, Área de Trabalho e pastas com nomes longos.

Use:

```powershell
cd C:\
mkdir dev
cd C:\dev
git clone URL_DO_REPOSITORIO tese
cd tese
```

## 4. Instalar camada Python

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements-dev.txt
python -m pip install -e .
python -m pytest -q
```

Se a ativação do ambiente virtual for bloqueada:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
.\.venv\Scripts\Activate.ps1
```
