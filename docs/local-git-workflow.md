# Fluxo Git Local do Projeto

Este projeto está sendo conduzido com fluxo local, sem Pull Request obrigatório.

## Regra

```text
branch → alterações → testes → commit → merge local → tag → push posterior
```

## Nova fase

```powershell
git checkout main
git status
git checkout -b tipo/nome-da-fase
```

Exemplos:

```powershell
git checkout -b docs/readme-thesis-summary
git checkout -b feat/case12-source-terms
git checkout -b chore/remove-generated-files
```

## Testes antes do commit

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
```

Quando houver visualização/dataset:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

## Commit

```powershell
git add -A
git commit -m "docs: add thesis README with equations and figures"
```

## Merge local

```powershell
git checkout main
git merge --no-ff docs/readme-thesis-summary -m "merge: add thesis README with equations and figures"
```

## Tag

```powershell
git tag -a v0.5.2 -m "v0.5.2: thesis README with equations and figures"
```

## Push posterior

```powershell
git push origin main
git push origin --tags
```

## Não usar agora

```powershell
git push --force
```

Sem necessidade. `--force` é a motosserra do Git: útil quando se sabe exatamente o que está cortando; perigoso quando só se queria podar uma roseira.
