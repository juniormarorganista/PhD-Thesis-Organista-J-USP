# Case 12 `tfwz` Fortran reference exporter

This directory contains a standalone Fortran program that exports reference
values for the MMS Case 12 vorticity source term:

```text
tfwz
```

The expression mirrors the original routine:

```text
src/fortran_original/4Order_LPOG_Euler/case_12.f::source_term_case_12(t)
```

The original Fortran thesis code is not modified.

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_tfwz\case12_tfwz_reference.f90 `
  -o validation\fortran_reference\case12_tfwz\case12_tfwz_reference.exe
```

## Generate reference CSV

Default `3 x 3` grid:

```powershell
.\validation\fortran_reference\case12_tfwz\case12_tfwz_reference.exe `
  validation\reference_data\case12_tfwz_3x3_default.csv 3 3
```

A denser `5 x 5` grid:

```powershell
.\validation\fortran_reference\case12_tfwz\case12_tfwz_reference.exe `
  validation\reference_data\case12_tfwz_5x5_default.csv 5 5
```

## Validate Python against reference

```powershell
python examples\05_validate_case12_tfwz_reference.py `
  validation\reference_data\case12_tfwz_3x3_default.csv
```

and:

```powershell
python examples\05_validate_case12_tfwz_reference.py `
  validation\reference_data\case12_tfwz_5x5_default.csv
```
