# GitHub workflow of the project

This project is conducted through pull requests on GitHub.

## Rule

```text
main → branch → commits → push → PR on GitHub → merge via GitHub → main
```

## New phase

```powershell
git checkout main
git pull --ff-only origin main
git checkout -b tipo/nome-da-fase
```

Examples:

```powershell
git checkout -b docs/readme-thesis-summary
git checkout -b feat/case12-source-terms
git checkout -b chore/remove-generated-files
```

## Tests before committing

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
```

When there is visualization/dataset work:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

## Commit and push

```powershell
git add <files>
git commit -m "docs: add thesis README with equations and figures"
git push -u origin docs/readme-thesis-summary
```

## Pull request

```powershell
gh pr create --base main --fill
gh pr checks --watch
gh pr merge --merge --delete-branch
```

## Not used

```powershell
git push --force
```

Not needed. `--force` is the chainsaw of Git: useful when you know exactly
what you are cutting; dangerous when you only wanted to prune a rose bush.
