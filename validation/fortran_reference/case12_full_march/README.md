# Case 12 full time march Fortran reference

This directory contains a standalone Fortran program that repeats the full
time-marching step from Etapa 10/11
(`validation/fortran_reference/case12_multistep`): explicit Euler step (Etapa
9) + boundary-condition re-imposition at `t+dt` + Poisson re-solve for `psi`
(Etapa 8) + recomputation of `u`/`v` from the updated `psi`, in a loop that
mirrors the *stopping criterion* of `solver_euler`
(`src/fortran_original/4Order_LPOG_Euler/euler_rk.f`, lines 326-417):

```text
do while (time < timef - dt)
  ... drv + Euler update + source_term + boundaries_condictions + loop(2) ...
  errowz  = max|wz_new  - wz_old | / max(max|wz_new |, 1)
  erroTxx = max|Txx_new - Txx_old| / max(max|Txx_new|, 1)
  erroTxy = max|Txy_new - Txy_old| / max(max|Txy_new|, 1)
  erroTyy = max|Tyy_new - Tyy_old| / max(max|Tyy_new|, 1)
  emv = max(errowz, erroTxx, erroTxy, erroTyy)
  if (emv < 1.0e-7 .and. time > 0.01) exit  ! converged to steady state
end do
```

The original Fortran thesis code is not modified.

## Scope

- Same second-order finite-difference operators, dense Gaussian-elimination
  Poisson re-solve and `wz`-boundary simplification as Etapas 10/11. See
  `validation/fortran_reference/case12_multistep/README.md` and
  `docs/case-12-multistep-marching-comparison.md`.
- Default parameters (`Re=1`, `Wi=1`, `beta=0.5`, `epsilon=0`, `xi=0`,
  `alphaG=0`, `dt=1.0e-5`, `at=0`) match
  `mms_viscoelastic.numerics.time_stepping.EulerStepParameters()`,
  `Case12Parameters()` and `Case12SourceTermParameters()`.
- Time-independence at `at = 0` (boundary fields and source terms reused from
  the `t=0` input CSV at every step) as in Etapas 10/11.

## Divergence detection

If any of `errowz`/`erroTxx`/`erroTxy`/`erroTyy`/`emv` becomes `NaN` (e.g.
because a field overflowed to `+/-Infinity`), the march stops immediately
with `diverged = .true.`. **This guard is important**: without it,
gfortran's `max(NaN, x) = x` semantics make `emv` collapse silently to `0.0`
once a `NaN` appears anywhere in the comparison, which would otherwise be
misreported as `converged = .true.`.

## Key finding (Etapa 12)

On a `9 x 9` grid, this march **does not converge** to `emv < 1e-7`. Instead
`emv` decreases for a short while, then grows (slowly at first, then
explosively) until the fields overflow at step 455064 (`time ~= 4.5506`),
which is correctly reported as `diverged = .true.`. See
[`docs/case-12-full-march-reference.md`](../../../docs/case-12-full-march-reference.md)
for the full discussion, data and the recommended "stable horizon" for Etapa
13.

## Generate the input CSV

```powershell
python validation\fortran_reference\case12_time_step\generate_input_csv.py `
  validation\reference_data\case12_time_step_input_9x9.csv 9 9

python validation\fortran_reference\case12_time_step\generate_input_csv.py `
  validation\reference_data\case12_time_step_input_17x17.csv 17 17
```

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_full_march\case12_full_march_reference.f90 `
  -o validation\fortran_reference\case12_full_march\case12_full_march_reference.exe
```

## Run

```text
case12_full_march_reference.exe <infile> <outfile> <nx> <ny> <dt> [histfile] [max_steps] [timef] [hist_start_step] [hist_stride]
```

- `histfile` (optional): if non-empty, writes a per-step CSV with columns
  `step,time,errowz,erroTxx,erroTxy,erroTyy,emv,maxwz,maxtxx,maxtxy,maxtyy,maxpsi`.
  The `max*` columns are the per-step maximum absolute field values, used to
  detect blow-up independently of the relative error metrics.
- `hist_start_step` / `hist_stride` (optional, default `1`/`1`): restrict the
  history file to rows with `t_step >= hist_start_step` and
  `(t_step - hist_start_step) mod hist_stride == 0`, keeping the file small
  for long runs.
- `max_steps` (default `200000`) is a safety cap independent of `solver_euler`'s
  own `time < timef - dt` bound (`timef = 50` would correspond to ~5e6 steps
  at `dt = 1e-5`).

Example (`9 x 9`, full history every 200 steps up to the divergence point):

```powershell
.\validation\fortran_reference\case12_full_march\case12_full_march_reference.exe `
  validation\reference_data\case12_time_step_input_9x9.csv `
  validation\reference_data\case12_full_march_output_9x9.csv `
  9 9 1.0e-5 `
  validation\reference_data\case12_full_march_history_9x9.csv `
  460000 50.0 1 200
```

Example (`17 x 17`, stable horizon up to `t = 1.5`):

```powershell
.\validation\fortran_reference\case12_full_march\case12_full_march_reference.exe `
  validation\reference_data\case12_time_step_input_17x17.csv `
  validation\reference_data\case12_full_march_output_17x17.csv `
  17 17 1.0e-5 `
  validation\reference_data\case12_full_march_history_17x17.csv `
  150000 50.0 1 500
```

## CSV schema

Input (`case12_time_step_input_<nx>x<ny>.csv`):

```text
x,y,u,v,wz,psi,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy
```

Output (`case12_full_march_output_<nx>x<ny>.csv`): field values at the last
completed step (last row before divergence, or `max_steps`/convergence if
reached):

```text
x,y,wz,txx,txy,tyy,psi,u,v
```

History (`case12_full_march_history_<nx>x<ny>.csv`):

```text
step,time,errowz,erroTxx,erroTxy,erroTyy,emv,maxwz,maxtxx,maxtxy,maxtyy,maxpsi
```
