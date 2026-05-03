# Release Notes — v0.6.7.1

## Summary

Fixes the validation layer to match the current `Grid2D.from_bounds` API.

## Fixed

- Removed unsupported `x0`, `x1`, `y0`, `y1` keyword arguments from validation
  grid reconstruction.
- Added explicit validation that Fortran reference CSV coordinates are uniform
  grids on `[0, 1]`.
- Added regression tests for reference validation grid reconstruction.

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```
