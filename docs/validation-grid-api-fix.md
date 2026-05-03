# Validation grid API fix

## Problem

The validation helper previously called:

```python
Grid2D.from_bounds(
    nx=nx,
    ny=ny,
    x0=...,
    x1=...,
    y0=...,
    y1=...,
)
```

However, the current `Grid2D.from_bounds` API accepts only:

```python
Grid2D.from_bounds(nx=nx, ny=ny)
```

This caused all Fortran-reference validation tests to fail with:

```text
TypeError: Grid2D.from_bounds() got an unexpected keyword argument 'x0'
```

## Fix

The helper `_grid_from_reference` now:

1. reads the unique `x` and `y` values from the reference CSV;
2. validates that they form uniform unit-square grids;
3. calls the supported API:

```python
Grid2D.from_bounds(nx=nx, ny=ny)
```

## Why this is safe

The current Fortran reference CSVs all use:

```text
x in [0, 1]
y in [0, 1]
```

with uniform spacing. If a future reference dataset uses other bounds, the
validator now fails explicitly instead of silently comparing against the wrong
grid.

## Validation

Run:

```powershell
python -m pytest -q
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```
