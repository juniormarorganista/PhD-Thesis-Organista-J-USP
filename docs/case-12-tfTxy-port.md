# MMS Case 12 — `tfTxy` Python port

This document records the third source-term port from the original Fortran file:

```text
src/fortran_original/4Order_LPOG_Euler/case_12.f
```

## Ported components

```text
tfwz
tfTxx
tfTxy
```

Still pending:

```text
tfTyy
```

## Python API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import evaluate_case_12_tfTxy

grid = Grid2D.from_bounds(nx=65, ny=65)
tfTxy = evaluate_case_12_tfTxy(grid)
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
terms.tfTxy
```

Pending:

```text
terms.tfTyy
```

## Example

```powershell
python examples\04_case12_tfwz_source_term.py
```

This exports the ported terms to:

```text
outputs/case12/source_terms/case12_ported_sources.csv
```

Current columns:

```text
x
y
tfwz
tfTxx
tfTxy
```

## Validation status

Current tests verify:

- registry status;
- shape convention;
- finite values;
- default reference value at the origin.

Still required:

- regression against Fortran-generated reference data for `tfTxy`;
- port of `tfTyy`;
- complete Case 12 source-term validation.
