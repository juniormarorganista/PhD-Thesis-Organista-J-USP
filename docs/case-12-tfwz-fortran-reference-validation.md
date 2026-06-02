# MMS Case 12 - Fortran reference validation for `tfwz`

This phase adds a small Fortran reference exporter and Python comparison tools
for the first ported source term:

```text
tfwz
```

## Why this matters

The Python function:

```python
evaluate_case_12_tfwz(...)
```

contains a long algebraic expression. Before porting `tfTxx`, `tfTxy` and
`tfTyy`, this expression must be checked against independent Fortran-generated
reference values.

## Reference files

The repository now includes:

```text
validation/reference_data/case12_tfwz_3x3_default.csv
validation/reference_data/case12_tfwz_5x5_default.csv
```

These files were generated using the standalone Fortran exporter in:

```text
validation/fortran_reference/case12_tfwz/case12_tfwz_reference.f90
```

## Validate

```powershell
python examples\05_validate_case12_tfwz_reference.py `
  validation\reference_data\case12_tfwz_3x3_default.csv
```

```powershell
python examples\05_validate_case12_tfwz_reference.py `
  validation\reference_data\case12_tfwz_5x5_default.csv
```

Or run the full test suite:

```powershell
python -m pytest -q
```

## Compile the Fortran exporter

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_tfwz\case12_tfwz_reference.f90 `
  -o validation\fortran_reference\case12_tfwz\case12_tfwz_reference.exe
```

## Regenerate reference files

```powershell
.\validation\fortran_reference\case12_tfwz\case12_tfwz_reference.exe `
  validation\reference_data\case12_tfwz_3x3_default.csv 3 3
```

```powershell
.\validation\fortran_reference\case12_tfwz\case12_tfwz_reference.exe `
  validation\reference_data\case12_tfwz_5x5_default.csv 5 5
```

## Current acceptance criterion

```text
max_abs_error < 1e-11
```

The comparison currently targets the default parameter set:

```text
Re = 1
beta = 0.5
at = 0
t = 0
all trigonometric tensor coefficients = 1
```

## Next step

After this validation passes locally, proceed to:

```text
v0.6.3 - port tfTxx
```
