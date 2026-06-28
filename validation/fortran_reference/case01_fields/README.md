# Case 01 fields and source terms Fortran reference

This directory contains a standalone Fortran program that exports reference
values for the MMS Case 01 exact fields and source terms:

```text
u, v, wz, psi, txx, txy, tyy, tfwz, tfTxx, tfTxy, tfTyy
```

The expressions mirror the original routines:

```text
src/fortran_original/4Order_LPOG_Euler/case_01.f::exact_solutions_case_01(t)
src/fortran_original/4Order_LPOG_Euler/case_01.f::source_term_case_01(t)
```

`boundaries_condictions_case_01` is not separately exported: it re-evaluates
the same exact-solution formulas restricted to the boundary nodes, with no
distinct expression of its own (see the module docstring in
`src/mms_viscoelastic/manufactured_solutions/case_01.py`).

Unlike Case 12, Case 01's polymeric stress field (`txx`, `txy`, `tyy`) is
spatially uniform: it reproduces the fixed targets `c1=1`, `c2=2`, `c3=3`
from `par.nn`, scaled by `(1 - beta)`.

Default parameters (`Rey=1`, `Wi=1`, `beta=0.5`, `epsilon=0`, `xi=0`,
`alphaG=0`, `at=0`) match `Case01Parameters`/`Case01SourceTermParameters` in
the Python port, not the transient `Rey`/`Wi`/`epsylon` values currently
sitting in `par.nn` (which reflect whichever MMS case `par.nn` was last
configured to run via `ty_mms`).

The original Fortran thesis code is not modified.

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case01_fields\case01_fields_reference.f90 `
  -o validation\fortran_reference\case01_fields\case01_fields_reference.exe
```

## Generate reference CSV

```powershell
.\validation\fortran_reference\case01_fields\case01_fields_reference.exe `
  validation\reference_data\case01_fields_5x5_default.csv 5 5

.\validation\fortran_reference\case01_fields\case01_fields_reference.exe `
  validation\reference_data\case01_fields_9x9_default.csv 9 9
```

## Validate Python against reference

```powershell
python -m pytest tests\test_case01_fields_fortran_reference.py -v
python examples\16_case01_fields_and_source_terms.py
```
