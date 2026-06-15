# Phase 4 - Case 12 visualization and Machine Learning datasets

This phase adds a practical layer on top of the Python port of manufactured case 12.
The goal is to turn the exact solution into two useful products:

1. visual-inspection figures for the thesis, presentation and README;
2. tabular datasets for future Machine Learning experiments.

## Files created

```text
src/mms_viscoelastic/visualization/__init__.py
src/mms_viscoelastic/visualization/fields.py
src/mms_viscoelastic/ml/datasets.py
examples/03_case12_visualization_and_ml_dataset.py
tests/test_case12_visualization.py
tests/test_case12_datasets.py
```

## Axis convention

The `plot_scalar_field` function explicitly validates that the array has shape:

```text
(ny, nx)
```

This avoids accidentally swapping axes between the Fortran convention and the
NumPy/Matplotlib convention. This point is critical for correctly comparing
numerical results with the thesis figures.

## How to generate figures and datasets

From the repository root, after installing the package in development mode:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

With a finer grid:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py --nx 129 --ny 129
```

With a custom output folder:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py --output-dir outputs\case12
```

## Expected outputs

```text
outputs/case12/figures/case12_panel.png
outputs/case12/figures/case12_u.png
outputs/case12/figures/case12_v.png
outputs/case12/figures/case12_psi.png
outputs/case12/figures/case12_wz.png
outputs/case12/figures/case12_txx.png
outputs/case12/figures/case12_txy.png
outputs/case12/figures/case12_tyy.png
outputs/case12/datasets/case12_full.csv
outputs/case12/datasets/case12_train.csv
outputs/case12/datasets/case12_validation.csv
outputs/case12/datasets/case12_metadata.txt
```

## Dataset columns

The target columns are:

```text
u, v, psi, wz, txx, txy, tyy
```

The input columns include:

```text
x, y, time, beta, x2, y2, xy, sin_pi_x, cos_pi_x, sin_pi_y, cos_pi_y
```

This format already allows training supervised-regression models for
manufactured fields and, in the future, comparing classical ML models, PINNs
and neural operators.

## Limitation

This phase does not yet use data produced by Fortran. It works with the exact
Python solution of `case_12`. The Fortran vs Python comparison should come in
its own phase, after the reference results have been exported in a
standardized format.
