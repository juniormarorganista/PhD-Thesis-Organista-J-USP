# Case 12 exact fields Fortran reference

This directory contains an isolated Fortran reference program for the MMS Case
12 exact manufactured fields:

```text
u, v, wz, psi, txx, txy, tyy
```

The program is a verbatim port of the formulas implemented in
`exact_solutions_case_12(t)` (`src/fortran_original/4Order_LPOG_Euler/case_12.f`),
which is the same subroutine used as the t=0 reference for
`evaluate_case_12_exact` in the Python implementation.

The program is intentionally separated from the original thesis Fortran code.
It does not modify the original files under:

```text
src/fortran_original/
```

## Compile

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_fields\case12_fields_reference.f90 `
  -o validation\fortran_reference\case12_fields\case12_fields_reference.exe
```

## Generate reference CSV

3 x 3:

```powershell
.\validation\fortran_reference\case12_fields\case12_fields_reference.exe `
  validation\reference_data\case12_fields_3x3_default.csv 3 3
```

5 x 5:

```powershell
.\validation\fortran_reference\case12_fields\case12_fields_reference.exe `
  validation\reference_data\case12_fields_5x5_default.csv 5 5
```

## CSV schema

```text
x,y,u,v,wz,psi,txx,txy,tyy
```

The grid uses the same convention as the Python validation layer:

```text
x in [0, 1]
y in [0, 1]
```

## Default parameters

The program uses the same defaults as `Case12Parameters()` in
`src/mms_viscoelastic/manufactured_solutions/case_12.py`:

- `betann = 0.5`
- `at = 0.0` (i.e. `t = 0`, the initial instant)
- All `a1..a3` and `b1..b3` trigonometric stress coefficients equal to `1.0`

## Scope

This reference covers only the **t=0 initial-instant** exact fields used by
`exact_solutions_case_12`. It does **not** reproduce the time-marching solver
loop (`source_term_case_12`, the discretized PDE update, etc.). A full
solver-level (multi time-step) comparison is left for a future phase; see
`docs/case-12-solver-comparison.md`.
