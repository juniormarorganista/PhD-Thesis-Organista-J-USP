# Case 12 multi-step time-marching Fortran reference

This directory contains a standalone Fortran program that exports reference
values after `n_steps` repetitions of the full time-marching step from Etapa
10 (`validation/fortran_reference/case12_time_step`): explicit Euler step
(Etapa 9) + boundary-condition re-imposition at `t+dt` + Poisson re-solve for
`psi` (Etapa 8) + recomputation of `u`/`v` from the updated `psi`.

Each iteration mirrors the sequence

```text
call source_term(time + dt)
call boundaries_condictions(time + dt)
call loop(2)
```

from `src/fortran_original/4Order_LPOG_Euler/euler_rk.f::solver_euler`. The
original Fortran thesis code is not modified.

## Scope

- Same second-order finite-difference operators, Poisson re-solve (dense
  Gaussian elimination vs `solve_poisson_dirichlet`'s sparse direct solve) and
  `wz`-boundary simplification as Etapa 10. See
  `validation/fortran_reference/case12_time_step/README.md` and
  `docs/case-12-time-marching-comparison.md`.
- Default parameters (`Re=1`, `Wi=1`, `beta=0.5`, `epsilon=0`, `xi=0`,
  `alphaG=0`, `dt=1.0e-5`, `at=0`) match
  `mms_viscoelastic.numerics.time_stepping.EulerStepParameters()`,
  `Case12Parameters()` and `Case12SourceTermParameters()`.

### Time-independence at `at = 0`

With the default manufactured-solution parameter `at = 0`, both the exact
fields `exact_solutions_case_12(t)` and the manufactured source terms
`tfwz`/`tfTxx`/`tfTxy`/`tfTyy` are independent of `t` (`exp(at*t) == 1` for
all `t`). This program therefore evaluates the boundary fields once and reuses
the `tfwz`/`tftxx`/`tftxy`/`tftyy` columns from the `t=0` input CSV at every
step, instead of regenerating per-step input CSVs.
`case12_multistep_march` (the Python counterpart) does not rely on this
property -- it re-evaluates the source terms and exact fields at the correct
`t` at every step, which is a no-op for `at = 0`. This is documented in
`docs/case-12-multistep-marching-comparison.md`.

## Generate the input CSV

The input CSV is the same one used by Etapa 10 (`t=0` exact fields and source
terms on a uniform grid):

```powershell
python validation\fortran_reference\case12_time_step\generate_input_csv.py `
  validation\reference_data\case12_time_step_input_9x9.csv 9 9
```

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_multistep\case12_multistep_reference.f90 `
  -o validation\fortran_reference\case12_multistep\case12_multistep_reference.exe
```

## Generate reference CSV

Default `9 x 9` grid, `dt = 1.0e-5`, `n_steps = 10`:

```powershell
.\validation\fortran_reference\case12_multistep\case12_multistep_reference.exe `
  validation\reference_data\case12_time_step_input_9x9.csv `
  validation\reference_data\case12_multistep_output_9x9_dt1e-05_n10.csv `
  9 9 1.0e-5 10
```

## CSV schema

Input (`case12_time_step_input_9x9.csv`):

```text
x,y,u,v,wz,psi,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy
```

Output (`case12_multistep_output_9x9_dt1e-05_n10.csv`):

```text
x,y,wz,txx,txy,tyy,psi,u,v
```

All fields in the output are the values **after** `n_steps` full
time-marching steps of size `dt` each.

## Validate Python against reference

```powershell
python examples\15_case12_multistep_comparison.py `
  validation\reference_data\case12_multistep_output_9x9_dt1e-05_n10.csv
```
