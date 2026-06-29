# Case 04 fields and source terms Fortran reference

This directory contains a standalone Fortran program that exports reference
values for the MMS Case 04 exact fields and source terms:

```text
u, v, wz, psi, txx, txy, tyy, tfwz, tfTxx, tfTxy, tfTyy
```

The expressions mirror the original routines:

```text
src/fortran_original/4Order_LPOG_Euler/case_04.f::exact_solutions_case_04(t)
src/fortran_original/4Order_LPOG_Euler/case_04.f::source_term_case_04(t)
```

`boundaries_condictions_case_04` is not separately exported, for the same
reason as Case 01 (see `validation/fortran_reference/case01_fields/README.md`).

`u`, `v`, `wz`, `psi` follow a Taylor-Green-vortex trigonometric ansatz
(`u = 8 cos(pi x) sin(pi y)`, etc.), distinct from cases 01/02/03 and case
12. `txx`, `txy`, `tyy` are spatially uniform constants, algebraically
identical to case 01's stress ansatz.

Default parameters (`Rey=1`, `Wi=1`, `beta=0.5`, `epsilon=0`, `xi=0`,
`alphaG=0`, `at=0`, `c1=1`, `c2=2`, `c3=3`) match
`Case04Parameters`/`Case04SourceTermParameters` in the Python port.

The original Fortran thesis code is not modified.

## Transcription method, and an error this caught

Generated with the same Fortran-to-Fortran mechanical script used for
Cases 01-03 (strips the Fortran77 continuation prefix, renames intrinsics,
preserves every parenthesis character-for-character).

Cross-validating this reference against the Python port
(`source_terms/case_04.py`) caught a real sign error: because Case 04's
source terms are short enough to read at a glance, the outer
`(-1 + betann)` factor multiplying `tfTxx`/`tfTxy`/`tfTyy` was retyped by
hand as `(1 - betann)` -- the negative of the correct value. The error
showed up as a constant ~100-200 absolute discrepancy at every grid
point (independent of `x`/`y`, the signature of a sign flip in a
position-independent term), and was fixed by restoring the original
`(-1 + betann)`. See `docs/case-04-port.md` for the full account.

## Compile

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case04_fields\case04_fields_reference.f90 `
  -o validation\fortran_reference\case04_fields\case04_fields_reference.exe
```

## Generate reference CSV

```powershell
.\validation\fortran_reference\case04_fields\case04_fields_reference.exe `
  validation\reference_data\case04_fields_5x5_default.csv 5 5

.\validation\fortran_reference\case04_fields\case04_fields_reference.exe `
  validation\reference_data\case04_fields_9x9_default.csv 9 9
```

## Validate Python against reference

```powershell
python -m pytest tests\test_case04_fields_fortran_reference.py -v
python examples\19_case04_fields_and_source_terms.py
```
