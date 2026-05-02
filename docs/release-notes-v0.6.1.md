# Release Notes — v0.6.1

## Summary

Ports the MMS Case 12 vorticity source term `tfwz` from Fortran to Python.

## Added

- `evaluate_case_12_tfwz`
- partial source-term support through `allow_partial=True`
- explicit registry of ported and pending Case 12 source-term components
- example script for exporting `tfwz` to CSV
- tests for shape, finite values, partial object behavior and a default reference value
- documentation for the `tfwz` port

## Still pending

- `tfTxx`
- `tfTxy`
- `tfTyy`
- direct regression against Fortran-generated reference output

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\04_case12_tfwz_source_term.py
```
