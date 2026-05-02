# Local Git Workflow

This project is currently using a local-maintainer workflow.

The preferred sequence is:

```text
branch → changes → tests → commit → local merge into main → tag → later push
```

## Start a new phase

```powershell
git checkout main
git pull --ff-only origin main
git checkout -b feat/name-of-phase
```

## Validate changes

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
```

Run relevant examples:

```powershell
python examples\01_smoke_python_package.py
python examples\02_case12_manufactured_solution.py
python examples\03_case12_visualization_and_ml_dataset.py
```

## Commit phase

```powershell
git status
git diff --stat
git add .
git commit -m "feat: short objective message"
```

## Merge locally

```powershell
git checkout main
git merge --no-ff feat/name-of-phase -m "merge: short objective message"
```

## Tag locally

```powershell
git tag -a v0.X.0 -m "v0.X.0: release summary"
```

## Synchronize later

```powershell
git push origin main
git push origin --tags
```

Optionally push feature branches:

```powershell
git push origin feat/name-of-phase
```

## Rules

Do not use force push unless intentionally rewriting history.

Avoid committing generated `outputs/` data. Use `assets/preview/` only for curated figures intended for README or documentation.
