# ML package structure

The `mms_viscoelastic.ml` package is split into two layers.

```text
src/mms_viscoelastic/ml/
├── __init__.py
├── features.py
└── datasets.py
```

## `features.py`

Low-level feature engineering and conversion utilities:

- `grid_coordinates`
- `field_mapping`
- `add_engineered_features`
- `field_to_feature_table`

This module does not import from `mms_viscoelastic.ml`, preventing circular imports.

## `datasets.py`

Higher-level dataset builders:

- `build_case12_dataset`
- `split_train_validation`
- `save_dataset_bundle`

This module imports `field_to_feature_table` from:

```python
from mms_viscoelastic.ml.features import field_to_feature_table
```

not from:

```python
from mms_viscoelastic.ml import field_to_feature_table
```

The second form may trigger circular imports while `mms_viscoelastic.ml.__init__` is still being initialized.

## Local validation

```powershell
python -m pip install -e .
python -m pytest -q
python examples\03_case12_visualization_and_ml_dataset.py
```

## Cleanup note

If an older file exists at:

```text
src/mms_viscoelastic/ml.py
```

remove it, because the ML layer is now a package directory:

```powershell
Remove-Item .\src\mms_viscoelastic\ml.py -Force -ErrorAction SilentlyContinue
```
