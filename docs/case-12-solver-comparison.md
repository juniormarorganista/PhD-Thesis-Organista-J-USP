# MMS Case 12 - solver comparison (Fortran x Python)

This document records the FASE VII comparison between the original Fortran
solver and the Python `mms_viscoelastic` implementation for Case 12.

## Scope

The planning document (`docs/planejamento-tecnico.md`) originally framed
FASE VII as a full solver-level comparison: run the Fortran time-marching
solver (`src/fortran_original/4Order_LPOG_Euler/`) and the equivalent Python
solver side by side, then compare the fields `u`, `v`, `wz`, `psi` and the
stress tensor `T` across the time evolution.

That comparison is **not implemented here**. The `mms_viscoelastic` package
currently provides:

- the exact manufactured fields (`evaluate_case_12_exact`);
- the manufactured source terms (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`);
- constitutive models (UCM, Oldroyd-B);
- a grid-convergence study of the finite-difference operators.

It does **not** provide a discretized time-marching solver (no Poisson
solver for `psi`, no explicit/implicit time integration of the vorticity and
stress transport equations). Porting that solver is a substantially larger
undertaking than any single phase completed so far and is left for a future
phase.

### What this phase covers instead

This phase implements the **t=0 initial-instant field comparison**: an
isolated Fortran reference program reproduces
`exact_solutions_case_12(t)` from
`src/fortran_original/4Order_LPOG_Euler/case_12.f` (lines 4-49) verbatim, and
the Python `evaluate_case_12_exact` fields are compared against it at
`t = 0` for the default Case 12 parameters (`betann = 0.5`, `at = 0.0`, all
trigonometric stress coefficients equal to `1.0`).

Both the Fortran subroutine and `evaluate_case_12_exact` evaluate the same
closed-form manufactured solution; this comparison therefore verifies that the
Python port of the exact fields is bit-for-bit equivalent (to floating-point
precision) with the Fortran formulas that seed every solver run, for the
seven fields used throughout the thesis:

```text
u, v, wz, psi, txx, txy, tyy
```

## Files

```text
validation/fortran_reference/case12_fields/case12_fields_reference.f90
validation/reference_data/case12_fields_3x3_default.csv
validation/reference_data/case12_fields_5x5_default.csv
src/mms_viscoelastic/validation/solver_comparison.py
examples/11_case12_solver_comparison.py
tests/test_case12_solver_comparison.py
```

## Run validation

```powershell
python examples\11_case12_solver_comparison.py
```

For the 5 x 5 reference:

```powershell
python examples\11_case12_solver_comparison.py `
  validation\reference_data\case12_fields_5x5_default.csv
```

Expected result: `passed=True` for all seven components.

## Rebuild reference data

If `gfortran` is installed:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_fields\case12_fields_reference.f90 `
  -o validation\fortran_reference\case12_fields\case12_fields_reference.exe

.\validation\fortran_reference\case12_fields\case12_fields_reference.exe `
  validation\reference_data\case12_fields_3x3_default.csv 3 3

.\validation\fortran_reference\case12_fields\case12_fields_reference.exe `
  validation\reference_data\case12_fields_5x5_default.csv 5 5
```

## Acceptance criterion

The Python implementation passes when, for every component in
`CASE12_FIELD_COMPONENTS = ("u", "v", "wz", "psi", "txx", "txy", "tyy")`:

```text
max_abs_error < 1e-10
max_rel_error < 1e-10
```

for both the 3x3 and 5x5 default reference datasets. Observed errors are at
floating-point round-off level (~1e-16 to 1e-21), as expected for two
implementations of the same closed-form expression.

## Deferred work

A full solver-level comparison (time-marching loop, Poisson solve for `psi`,
vorticity and stress transport, multiple time steps) requires porting the
discretized solver from `src/fortran_original/4Order_LPOG_Euler/` to Python.
This is out of scope for the current phase and is tracked as future work in
`docs/planejamento-tecnico.md`.
