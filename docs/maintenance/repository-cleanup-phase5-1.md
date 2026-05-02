# Repository Cleanup Report — Phase 5.1

- Generated at: `2026-05-02 11:01:17`
- Branch: `chore/cleanup-generated-files`
- Mode: `apply`

## Rule

This cleanup phase preserves original scientific sources:

- `src/fortran_original/`
- `src/mathematica_original/`
- `thesis/original/`
- `papers/original/`

No Fortran, LaTeX, Mathematica or original thesis/paper file should be removed by this script.

## Planned/applied actions

- No cleanup actions detected.

## Git status at report time

```text
M .gitignore
 D README.before-phase2.md
 D README.before-phase5.md
 D criar_readme_visual_fase5.py
?? docs/archive/
?? docs/maintenance/
?? scripts/cleanup_repository_phase5_1.py
?? scripts/generators/
```

## Recommended validation

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\03_case12_visualization_and_ml_dataset.py
git status
git diff --stat
```
