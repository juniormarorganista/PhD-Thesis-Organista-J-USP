# Case 07 fields and source terms Fortran reference

This directory contains a standalone Fortran program that exports reference
values for the MMS Case 07 exact fields and source terms:

```text
u, v, wz, psi, txx, txy, tyy, tfwz, tfTxx, tfTxy, tfTyy
```

The expressions mirror the original routines:

```text
src/fortran_original/4Order_LPOG_Euler/case_07.f::exact_solutions_case_07(t)
src/fortran_original/4Order_LPOG_Euler/case_07.f::source_term_case_07(t)
```

`boundaries_condictions_case_07` is not separately exported, for the same
reason as Case 01 (see `validation/fortran_reference/case01_fields/README.md`).

Unlike every other ported case, `u`, `v`, `wz`, `psi` decay in time
through a hyperbolic tangent (`-1 + tanh(4 - 8*at*t)`) rather than
`exp(-at*t)`. `txx`, `txy`, `tyy` are spatially uniform constants,
algebraically identical to Case 01's stress ansatz.

Default parameters (`Rey=1`, `Wi=1`, `beta=0.5`, `epsilon=0`, `xi=0`,
`alphaG=0`, `at=0`, `c1=1`, `c2=2`, `c3=3`) match
`Case07Parameters`/`Case07SourceTermParameters` in the Python port.

The original Fortran thesis code is not modified.

## New intrinsics and literal syntax in this case

`case_07.f` introduces two things the mechanical transcription scripts
had not encountered in Cases 01-06:

- `Tanh`/`Cosh` (hyperbolic functions) -> renamed to `np.tanh`/`np.cosh`
  in the Python script (Fortran's generic `tanh`/`cosh` need no rename).
- Fortran double-precision literal exponents, e.g. `8.0d0`, `1.d0` ->
  the Python script now converts the `d`/`D` exponent marker to `e` (so
  `8.0d0` becomes the valid Python literal `8.0e0`), only when it
  immediately follows a digit or decimal point (so the `d` in `dexp`,
  which starts at a word boundary, is never touched).

## A floating-point evaluation-order artifact (not a transcription bug)

`tfTyy`'s mechanically-extracted expression, when compiled as a single
deeply nested Fortran statement, evaluated to a value ~4.5e-8 away from
the value obtained either in Python or in a Fortran re-implementation
that assigns the same sub-expressions to named intermediate variables
(`term_a`, `term_b`, `term_c`, `inner`, `numer`, `denom`) before
combining them. This was confirmed to be a floating-point
evaluation-order/cancellation artifact specific to gfortran's compilation
of that one statement -- not a transcription error -- by:

1. Verifying the monolithic Fortran statement is byte-identical (modulo
   whitespace) to the mechanically extracted text from `case_07.f`.
2. Verifying the same text, evaluated in Python/NumPy, agrees with the
   decomposed Fortran re-implementation to ~1e-16 (ordinary rounding).
3. Reproducing the ~4.5e-8 discrepancy with both `-O0` and `-O2`, and
   with `-ffp-contract=off`, ruling out compiler optimization and FMA
   contraction as the cause.

`tftyy_case07` in this file therefore uses the decomposed,
named-intermediate-variable form rather than the literal one-statement
mechanical output. This is the same decomposition style already used
throughout the Python port (`term_1`, `term_2`, ...), applied here for
numerical robustness rather than transcription safety.

## Compile

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case07_fields\case07_fields_reference.f90 `
  -o validation\fortran_reference\case07_fields\case07_fields_reference.exe
```

## Generate reference CSV

```powershell
.\validation\fortran_reference\case07_fields\case07_fields_reference.exe `
  validation\reference_data\case07_fields_5x5_default.csv 5 5

.\validation\fortran_reference\case07_fields\case07_fields_reference.exe `
  validation\reference_data\case07_fields_9x9_default.csv 9 9
```

## Validate Python against reference

```powershell
python -m pytest tests\test_case07_fields_fortran_reference.py -v
python examples\22_case07_fields_and_source_terms.py
```
