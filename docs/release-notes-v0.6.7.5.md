# Release Notes — v0.6.7.5

## Summary

Adds thesis figures directly to the README using curated PNG assets converted from the thesis/slides EPS figures.

## Added

- `assets/thesis/case12/` with selected thesis figures;
- exact maps for Case 12 fields;
- exact surfaces for representative fields;
- convergence curves from thesis material;
- numerical slice figures;
- `docs/thesis-figures-used-in-readme.md` with source mapping;
- expanded README section: `Thesis figures in this README`.

## Policy

Do not commit every thesis figure. Keep `assets/thesis/` curated and README-facing.

## Validation

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
```
