# Release Notes — v0.6.2

## Summary

Adds Fortran-reference validation for the ported MMS Case 12 `tfwz` source term.

## Added

- standalone Fortran reference exporter for `tfwz`;
- reference CSV files for `3 x 3` and `5 x 5` grids;
- Python validation utility for comparing `tfwz` against reference CSVs;
- example CLI validation script;
- regression tests against Fortran-generated reference data;
- documentation for regenerating and validating the reference files.

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\05_validate_case12_tfwz_reference.py
```

## Next

`v0.6.3` should port `tfTxx` only after `tfwz` reference validation passes.
