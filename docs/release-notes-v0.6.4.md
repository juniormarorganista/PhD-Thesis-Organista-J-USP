# Release Notes — v0.6.4

## Summary

Adds Fortran reference validation for the MMS Case 12 tensor source term `tfTxx`.

## Added

- isolated Fortran reference program for `tfTxx`;
- 3 x 3 and 5 x 5 reference CSV files;
- Python validation utilities for source-term CSV comparison;
- example validation command;
- pytest regression tests for `tfTxx`;
- documentation for rebuilding reference data.

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\06_validate_case12_tfTxx_reference.py
```

## Next release

`v0.6.5` should port `tfTxy`.
