# Release Notes — v0.6.7.2

## Summary

Restores grid-shape metadata in source-term reference comparison objects.

## Fixed

- Added `nx` and `ny` to `SourceTermReferenceComparison`.
- Updated reference-grid reconstruction to return `(grid, nx, ny)`.
- Restored compatibility with existing `tfwz` tests and validation example.

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```
