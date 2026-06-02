# MMS Case 12 - `tfwz` Python port

This document records the first algebraic source-term port from the original
Fortran file:

```text
src/fortran_original/4Order_LPOG_Euler/case_12.f
```

## Ported component

```text
tfwz
```

The remaining tensor source terms are still pending:

```text
tfTxx
tfTxy
tfTyy
```

## Python API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import evaluate_case_12_tfwz

grid = Grid2D.from_bounds(nx=65, ny=65)
tfwz = evaluate_case_12_tfwz(grid)
```

The returned array uses NumPy convention:

```text
(ny, nx)
```

## Partial source-term object

A partial source-term object can be requested explicitly:

```python
from mms_viscoelastic.source_terms import evaluate_case_12_source_terms

terms = evaluate_case_12_source_terms(grid, allow_partial=True)
```

Then:

```text
terms.tfwz   -> NumPy array
terms.tfTxx  -> None
terms.tfTxy  -> None
terms.tfTyy  -> None
```

The full evaluator still raises `NotImplementedError` by default because the
complete source-term system has not been ported yet.

## Example

```powershell
python examples\04_case12_tfwz_source_term.py
```

This creates:

```text
outputs/case12/source_terms/case12_tfwz.csv
```

## Validation status

Current tests verify:

- import stability;
- source-term registry;
- shape convention;
- finite values;
- a default value at the origin.

Still required:

- regression against Fortran-generated reference data;
- tolerance-based comparison for selected grids;
- port of `tfTxx`;
- port of `tfTxy`;
- port of `tfTyy`.

## Next step

The next recommended phase is:

```text
v0.6.2 - port tfTxx
```

However, before porting the tensor source terms, it is better to create a small
Fortran reference-output workflow for `tfwz`, so the Python expression can be
validated against the original code.
