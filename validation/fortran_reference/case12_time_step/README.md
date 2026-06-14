# Case 12 time-marching step Fortran reference

This directory contains a standalone Fortran program that exports reference
values for one full time-marching step of the MMS Case 12 ``ty_sim == 1``
formulation: the explicit Euler step (Etapa 9) followed by boundary-condition
re-imposition at ``t+dt``, the Poisson re-solve for ``psi`` (Etapa 8) and the
recomputation of ``u``/``v`` from the updated ``psi``.

This mirrors the sequence

```text
call source_term(time + dt)
call boundaries_condictions(time + dt)
call loop(2)
```

from ``src/fortran_original/4Order_LPOG_Euler/euler_rk.f::solver_euler``. The
original Fortran thesis code is not modified.

## Scope

- Derivatives use the same second-order finite-difference operators as
  `mms_viscoelastic.numerics.derivatives`, **not** the 4th-order LPOG compact
  operators (`src/fortran_original/4Order_LPOG_Euler/derivs.f`) nor the Lele
  compact filter (`filter_trid_x`/`filter_trid_y` in `loop.f`). Porting those
  remains deferred.
- The Poisson re-solve uses a dense Gaussian elimination with partial
  pivoting over the same 5-point second-order Laplacian system assembled by
  `mms_viscoelastic.numerics.poisson.solve_poisson_dirichlet` (which uses a
  sparse direct solve). Both are exact solves of the same discrete linear
  system, so they agree to round-off.
- **`wz` boundary condition**: re-imposed from the exact manufactured solution
  at `t+dt`, *not* from the 4th-order one-sided stencil in
  `boundarie_vorticity_cavity_flow.f`. This is a deliberate deviation from the
  Fortran reference; see `docs/case-12-time-marching-comparison.md`.
- Default parameters (`Re=1`, `Wi=1`, `beta=0.5`, `epsilon=0`, `xi=0`,
  `alphaG=0`, `dt=1.0e-5`) match
  `mms_viscoelastic.numerics.time_stepping.EulerStepParameters()`.

## Generate the input CSV

The input CSV provides the `t=0` exact fields (`u`, `v`, `wz`, `psi`, `txx`,
`txy`, `tyy`) and source terms (`tfwz`, `tftxx`, `tftxy`, `tftyy`) from the
already-validated Python implementations:

```powershell
python validation\fortran_reference\case12_time_step\generate_input_csv.py `
  validation\reference_data\case12_time_step_input_9x9.csv 9 9
```

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_time_step\case12_time_step_reference.f90 `
  -o validation\fortran_reference\case12_time_step\case12_time_step_reference.exe
```

## Generate reference CSV

Default `9 x 9` grid with `dt = 1.0e-5`:

```powershell
.\validation\fortran_reference\case12_time_step\case12_time_step_reference.exe `
  validation\reference_data\case12_time_step_input_9x9.csv `
  validation\reference_data\case12_time_step_output_9x9_dt1e-05.csv `
  9 9 1.0e-5
```

## CSV schema

Input (`case12_time_step_input_9x9.csv`):

```text
x,y,u,v,wz,psi,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy
```

Output (`case12_time_step_output_9x9_dt1e-05.csv`):

```text
x,y,wz,txx,txy,tyy,psi,u,v
```

All fields in the output are the values **after** one full time-marching step
of size `dt` (Euler update + boundary re-imposition + Poisson re-solve +
velocity recomputation).

## Validate Python against reference

```powershell
python examples\14_case12_time_step_comparison.py `
  validation\reference_data\case12_time_step_output_9x9_dt1e-05.csv
```
