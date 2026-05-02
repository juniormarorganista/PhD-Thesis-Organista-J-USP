# Release Notes — v0.5.4

## Summary

Fixes compatibility issues in the Case 12 dataset and visualization layer.

## Fixed

- `Grid2D` is now instantiated with `Grid2D.from_bounds(nx=..., ny=...)`.
- `field_to_feature_table` now preserves the original conservative default:
  coordinates plus physical fields only.
- Engineered ML features are added only when explicitly requested.
- Visualization now exports backward-compatible names:
  - `plot_field`
  - `plot_field_panel`
- The Case 12 visualization example now closes Matplotlib figures after saving.

## Validation

```powershell
python -m pip install -e .
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\03_case12_visualization_and_ml_dataset.py
```
