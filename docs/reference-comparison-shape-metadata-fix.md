# Reference comparison shape metadata fix

## Problem

The reference comparison object was missing:

```text
nx
ny
```

Older tests and `examples/05_validate_case12_tfwz_reference.py` expect these
attributes:

```python
comparison.nx
comparison.ny
```

The numerical comparison itself was already passing, but the object schema was
incomplete.

## Fix

`SourceTermReferenceComparison` now includes:

```python
nx: int
ny: int
```

The helper `_grid_from_reference` now returns:

```python
(grid, nx, ny)
```

and `compare_case12_reference_csv` stores this metadata in the comparison
object.

## Validation

```powershell
python -m pytest -q
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```
