# Lele compact filter Fortran reference

This directory contains a standalone Fortran program that exports reference
values for the Lele-type compact spatial filter:

```text
lele_filter_x / lele_filter_y
```

The algorithm mirrors the original routines:

```text
src/fortran_original/4Order_LPOG_Euler/loop.f     (rhs_tridfx_wz, rhs_tridfy_wz, lhs_tridfx, lhs_tridfy)
src/fortran_original/4Order_LPOG_Euler/solvers.f  (tridx, tridy)
src/fortran_original/4Order_LPOG_Euler/par.nn     (alphaf = 0.48, af, bf, cf, df)
```

applied to `f(x, y) = sin(3*pi*x)*cos(2*pi*y) + 0.1*sin(20*pi*x)`, chosen to
have a high-frequency component so that the filter has a non-trivial effect
to check. The original Fortran thesis code is not modified.

Note: `rhs_tridfx_T` (used for `Txx`/`Txy`/`Tyy`) applies the identical
weights as `rhs_tridfx_wz`; this reference and the Python port
(`mms_viscoelastic.numerics.filters`) therefore implement the formula once,
generically, rather than once per field.

The doxygen comment above `rhs_tridfx_wz` in the original Fortran states
the boundary formula for `f_1` with both `-4f_4` and `-f_5` (i.e. all minus
signs); the executable code instead computes `+4f_4-f_5`. This reference
follows the executable code, not the comment.

## Compile

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case_lele_filter\case_lele_filter_reference.f90 `
  -o validation\fortran_reference\case_lele_filter\case_lele_filter_reference.exe
```

## Generate reference CSV

```powershell
.\validation\fortran_reference\case_lele_filter\case_lele_filter_reference.exe `
  validation\reference_data\case_lele_filter_13x13_default.csv 13 13

.\validation\fortran_reference\case_lele_filter\case_lele_filter_reference.exe `
  validation\reference_data\case_lele_filter_21x21_default.csv 21 21
```

## Validate Python against reference

```powershell
python -m pytest tests\test_lele_filter_fortran_reference.py -v
```
