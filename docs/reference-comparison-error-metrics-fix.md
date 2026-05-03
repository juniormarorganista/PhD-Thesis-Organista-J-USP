# Reference comparison error metrics fix

## Problem

After restoring `nx` and `ny`, the old `tfwz` validation example still failed
because it also prints:

```python
comparison.mean_abs_error
```

The comparison object had maximum errors but not mean errors.

## Fix

`SourceTermReferenceComparison` now includes:

```python
mean_abs_error: float
mean_rel_error: float
```

and `compare_case12_reference_csv` computes:

```python
mean_abs_error = mean(abs(predicted - reference))
mean_rel_error = mean(relative_error)
```

## Validation

```powershell
python -m pytest -q
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```
