# Release Notes — v0.6.3

## Summary

Ports the MMS Case 12 tensor source term `tfTxx` from Fortran to Python.

## Added

- `evaluate_case_12_tfTxx`
- updated partial source-term object containing `tfwz` and `tfTxx`
- updated source-term registry
- tests for `tfTxx` shape, finite values and a default reference value
- documentation for the `tfTxx` port

## Still pending

- `tfTxy`
- `tfTyy`
- direct regression against Fortran-generated reference output for `tfTxx`

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\04_case12_tfwz_source_term.py
```
