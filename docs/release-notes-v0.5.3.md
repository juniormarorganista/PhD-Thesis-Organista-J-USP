# Release Notes — v0.5.3

## Summary

Fixes the ML package architecture by separating low-level feature extraction from dataset builders.

## Fixed

- Circular import in `mms_viscoelastic.ml`.
- Dataset generation now imports `field_to_feature_table` from `mms_viscoelastic.ml.features`.
- Added regression tests for ML package imports.

## Added

- `src/mms_viscoelastic/ml/features.py`
- updated `src/mms_viscoelastic/ml/datasets.py`
- updated `src/mms_viscoelastic/ml/__init__.py`
- `tests/test_ml_imports.py`
- updated dataset tests
- documentation for the ML package structure

## Validation

```powershell
python -m pip install -e .
python -m pytest -q
python examples\03_case12_visualization_and_ml_dataset.py
```
