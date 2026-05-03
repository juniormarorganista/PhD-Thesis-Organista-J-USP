# MMS Case 12 — `tfTxx` Python port

This document records the second source-term port from the original Fortran
file:

```text
src/fortran_original/4Order_LPOG_Euler/case_12.f
```

## Ported components

```text
tfwz
tfTxx
```

Still pending:

```text
tfTxy
tfTyy
```

## Python API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import evaluate_case_12_tfTxx

grid = Grid2D.from_bounds(nx=65, ny=65)
tfTxx = evaluate_case_12_tfTxx(grid)
```

The returned array uses NumPy convention:

```text
(ny, nx)
```

## Partial source-term object

```python
from mms_viscoelastic.source_terms import evaluate_case_12_source_terms

terms = evaluate_case_12_source_terms(grid, allow_partial=True)
```

Available:

```text
terms.tfwz
terms.tfTxx
```

Pending:

```text
terms.tfTxy
terms.tfTyy
```

## Example

```powershell
python examples\04_case12_tfwz_source_term.py
```

This now exports both `tfwz` and `tfTxx` to:

```text
outputs/case12/source_terms/case12_ported_sources.csv
```

## Validation status

Current tests verify:

- registry status;
- shape convention;
- finite values;
- default reference value at the origin.

Still required:

- regression against Fortran-generated reference data;
- tolerance-based comparison for selected grids;
- port of `tfTxy`;
- port of `tfTyy`.

## Next recommended step

Before porting `tfTxy`, generate a Fortran reference CSV for `tfTxx`, as was done
for `tfwz`.
