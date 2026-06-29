# Case 09 fields and source terms Fortran reference

This directory contains a standalone Fortran program that exports reference
values for the MMS Case 09 exact fields and source terms:

```text
u, v, wz, psi, txx, txy, tyy, tfwz, tfTxx, tfTxy, tfTyy
```

The expressions mirror the original routines:

```text
src/fortran_original/4Order_LPOG_Euler/case_09.f::exact_solutions_case_09(t)
src/fortran_original/4Order_LPOG_Euler/case_09.f::source_term_case_09(t)
```

`boundaries_condictions_case_09` is not separately exported, for the same
reason as Case 01 (see `validation/fortran_reference/case01_fields/README.md`).

`u`, `v`, `wz`, `psi` are algebraically identical to Case 08. The
polymeric stress field uses the same trigonometric ansatz as Case
03/06/12.

Default parameters (`Rey=1`, `Wi=1`, `beta=0.5`, `epsilon=0`, `xi=0`,
`alphaG=0`, `at=0`, all `a`/`b` coefficients `=1`) match
`Case09Parameters`/`Case09SourceTermParameters` in the Python port.

The original Fortran thesis code is not modified.

## The largest ported case

`case_09.f` is the largest and densest case in the entire project (790
lines). Combining Case 08's explicit nonlinear convective Jacobian with
Case 06's trigonometric stress ansatz makes `tfTxy` the densest single
expression ported in the whole project: 296 matched parenthesis pairs.

As with Case 08, no part of either this reference or
`source_terms/case_09.py` was retyped by hand. Both were assembled
entirely from the mechanical scripts' output (the Python module via an
assembly script, then `ruff format` for whitespace only). Both
transcriptions agreed at machine precision on the first attempt, despite
the unprecedented scale.

## Compile

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case09_fields\case09_fields_reference.f90 `
  -o validation\fortran_reference\case09_fields\case09_fields_reference.exe
```

## Generate reference CSV

```powershell
.\validation\fortran_reference\case09_fields\case09_fields_reference.exe `
  validation\reference_data\case09_fields_5x5_default.csv 5 5

.\validation\fortran_reference\case09_fields\case09_fields_reference.exe `
  validation\reference_data\case09_fields_9x9_default.csv 9 9
```

## Validate Python against reference

```powershell
python -m pytest tests\test_case09_fields_fortran_reference.py -v
python examples\24_case09_fields_and_source_terms.py
```
