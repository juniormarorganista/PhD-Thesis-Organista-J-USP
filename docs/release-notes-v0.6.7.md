# Release Notes — v0.6.7

## Summary

Ports the MMS Case 12 tensor source term `tfTyy` from Fortran to Python.

## Added

- `evaluate_case_12_tfTyy`
- complete `Case12SourceTerms` object containing:
  - `tfwz`
  - `tfTxx`
  - `tfTxy`
  - `tfTyy`
- updated source-term registry
- tests for `tfTyy` shape, finite values and a default reference value
- documentation for the `tfTyy` port

## Status

All four Case 12 source terms are now implemented in Python.

Still pending:

- direct regression against Fortran-generated reference output for `tfTyy`;
- complete Case 12 validation report.

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\04_case12_tfwz_source_term.py
```
