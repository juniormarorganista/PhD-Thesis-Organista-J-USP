# Case 12 - Preliminary Results and ML Dataset

This document summarizes the current Python-side status of MMS Case 12.

## Status

Implemented:

- exact manufactured fields;
- structured 2D grid evaluation;
- visualization utilities;
- ML-ready tabular dataset generation;
- train/validation split;
- automated tests.

Implementado desde v0.6.7 / PR #5:

- port completo dos quatro termos fonte do Fortran (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`);
- validação de cada termo contra CSVs gerados pelo Fortran original (max\_abs\_error < 1e-10).

Ainda pendente:

- equivalência a nível de solver Fortran × Python para o Case 12;
- verificação de ordem de convergência do solver Python.

## Generated fields

```text
u
v
psi
wz
txx
txy
tyy
```

## Dataset design

The dataset currently uses pointwise information from a structured grid.

Inputs/features:

```text
x, y, time, beta,
x2, y2, xy,
sin_pi_x, cos_pi_x,
sin_pi_y, cos_pi_y
```

Targets:

```text
u, v, psi, wz, txx, txy, tyy
```

## Preview assets copied to repository

- `assets/preview/case12_panel.png`
- `assets/preview/case12_u.png`
- `assets/preview/case12_v.png`
- `assets/preview/case12_psi.png`
- `assets/preview/case12_wz.png`
- `assets/preview/case12_txx.png`
- `assets/preview/case12_txy.png`
- `assets/preview/case12_tyy.png`

## Generate outputs

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

Fine grid:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py --nx 129 --ny 129
```

## Copy figures to README assets

After generating the outputs, run:

```powershell
python scripts\criar_readme_visual_fase5.py
```

## Scientific use

This case is useful for:

- checking shape conventions;
- checking x/y plotting orientation;
- creating ML-ready supervised datasets;
- preparing future Fortran/Python equivalence tests;
- producing figures for README and thesis-presentation material.
