# Case 12 explicit Euler step Fortran reference

This directory contains a standalone Fortran program that exports reference
values for one explicit Euler step of the MMS Case 12 ``ty_sim == 1``
formulation (``wz``, ``Txx``, ``Txy``, ``Tyy``).

The right-hand sides mirror the original routine:

```text
src/fortran_original/4Order_LPOG_Euler/drvs.f::drv
```

and the Euler update performed by:

```text
src/fortran_original/4Order_LPOG_Euler/euler_rk.f::solver_euler
```

The original Fortran thesis code is not modified.

## Scope

- Derivatives use the same second-order finite-difference operators as
  `mms_viscoelastic.numerics.derivatives` (central differences with one-sided
  edges for first derivatives, 3-point central with zero boundary for second
  derivatives), **not** the 4th-order LPOG compact operators
  (`src/fortran_original/4Order_LPOG_Euler/derivs.f`) nor the Lele compact
  filter (`filter_trid_x`/`filter_trid_y` in `loop.f`). Porting those is
  deferred.
- Only interior grid points are advanced by one Euler step of size `dt`;
  boundary values are copied unchanged. Re-imposing boundary conditions and
  recomputing `psi`/`ux`/`uy` from the updated `wz`
  (`boundaries_condictions` + `loop` in the Fortran code) is left to the
  multi-step time-marching solver.
- Default parameters (`Re=1`, `Wi=1`, `beta=0.5`, `epsilon=0`, `xi=0`,
  `alphaG=0`, `dt=1.0e-5`) match
  `mms_viscoelastic.numerics.time_stepping.EulerStepParameters()`. See
  `docs/case-12-euler-step-comparison.md` for the full scope statement.

## Generate the input CSV

The input CSV provides the `t=0` exact fields (`u`, `v`, `wz`, `txx`, `txy`,
`tyy`) and source terms (`tfwz`, `tftxx`, `tftxy`, `tftyy`) from the
already-validated Python implementations, avoiding a second transcription of
the large `tfTxx`/`tfTxy`/`tfTyy` symbolic expressions into Fortran:

```powershell
python validation\fortran_reference\case12_euler_step\generate_input_csv.py `
  validation\reference_data\case12_euler_step_input_9x9.csv 9 9
```

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_euler_step\case12_euler_step_reference.f90 `
  -o validation\fortran_reference\case12_euler_step\case12_euler_step_reference.exe
```

## Generate reference CSV

Default `9 x 9` grid with `dt = 1.0e-5`:

```powershell
.\validation\fortran_reference\case12_euler_step\case12_euler_step_reference.exe `
  validation\reference_data\case12_euler_step_input_9x9.csv `
  validation\reference_data\case12_euler_step_output_9x9_dt1e-05.csv `
  9 9 1.0e-5
```

## CSV schema

Input (`case12_euler_step_input_9x9.csv`):

```text
x,y,u,v,wz,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy
```

Output (`case12_euler_step_output_9x9_dt1e-05.csv`):

```text
x,y,wz,txx,txy,tyy
```

`wz`, `txx`, `txy` and `tyy` in the output are the fields **after** one
explicit Euler step of size `dt` (interior points only; boundary values equal
the input).

## Validate Python against reference

```powershell
python examples\13_case12_euler_step_comparison.py `
  validation\reference_data\case12_euler_step_output_9x9_dt1e-05.csv
```
