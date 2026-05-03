# Release Notes — v0.6.5

## Summary

Ports the MMS Case 12 tensor source term `tfTxy` from Fortran to Python.

## Added

- `evaluate_case_12_tfTxy`
- updated partial source-term object containing `tfwz`, `tfTxx` and `tfTxy`
- updated source-term registry
- tests for `tfTxy` shape, finite values and a default reference value
- documentation for the `tfTxy` port

## Fixed

- Corrected the default test reference for `tfTxx` at the origin to `2.0`.

## Still pending

- `tfTyy`
- direct regression against Fortran-generated reference output for `tfTxy`

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\04_case12_tfwz_source_term.py
```
