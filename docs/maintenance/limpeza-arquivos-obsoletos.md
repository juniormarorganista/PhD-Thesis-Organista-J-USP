# Limpeza de arquivos obsoletos

Esta fase remove os scripts geradores/corretivos e backups antigos. A partir de agora, as entregas serão feitas por arquivos finais dentro de ZIP, sem scripts que geram outros arquivos.

## Arquivos a remover

Remover da raiz:

```text
criar_*.py
corrigir_*.py
reorganizar_*.ps1
cleanup_repository_phase*.py
README.before-*.md
README.backup*.md
README.old*.md
Texto colado.txt
```

Remover diretórios de scripts geradores/corretivos, se existirem:

```text
scripts/generators/
scripts/fixes/
```

## Comandos PowerShell

Execute na raiz do projeto:

```powershell
cd C:\dev\Thesis_USP_Codes_Organista_J

Remove-Item .\criar_*.py -Force -ErrorAction SilentlyContinue
Remove-Item .\corrigir_*.py -Force -ErrorAction SilentlyContinue
Remove-Item .\reorganizar_*.ps1 -Force -ErrorAction SilentlyContinue
Remove-Item .\cleanup_repository_phase*.py -Force -ErrorAction SilentlyContinue

Remove-Item .\README.before-*.md -Force -ErrorAction SilentlyContinue
Remove-Item .\README.backup*.md -Force -ErrorAction SilentlyContinue
Remove-Item .\README.old*.md -Force -ErrorAction SilentlyContinue
Remove-Item ".\Texto colado.txt" -Force -ErrorAction SilentlyContinue

Remove-Item .\scripts\generators -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item .\scripts\fixes -Recurse -Force -ErrorAction SilentlyContinue
```

## Não remover

Não remover:

```text
src/fortran_original/
src/mathematica_original/
thesis/original/
papers/original/
src/mms_viscoelastic/
tests/
examples/
docs/
assets/
```

## Depois da limpeza

```powershell
git status
git diff --stat

python -m compileall src\mms_viscoelastic
python -m pytest -q
```
