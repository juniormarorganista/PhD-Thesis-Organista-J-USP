# Release Notes — v0.6.7.3

## Summary

Restores mean error metrics in source-term reference comparison objects.

## Fixed

- Added `mean_abs_error` to `SourceTermReferenceComparison`.
- Added `mean_rel_error` to `SourceTermReferenceComparison`.
- Restored compatibility with the existing `tfwz` validation example.

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```
