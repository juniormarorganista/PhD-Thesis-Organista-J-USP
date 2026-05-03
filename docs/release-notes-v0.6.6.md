# Release Notes — v0.6.6

## Summary

Adds Fortran reference validation for the MMS Case 12 tensor source term `tfTxy`.

## Added

- isolated Fortran reference program for `tfTxy`;
- 3 x 3 and 5 x 5 reference CSV files;
- Python validation utility for `tfTxy`;
- example validation command;
- pytest regression tests for `tfTxy`;
- documentation for rebuilding reference data.

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\07_validate_case12_tfTxy_reference.py
```
