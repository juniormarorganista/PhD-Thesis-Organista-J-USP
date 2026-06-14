# MMS Case 12 - Explicit Euler step for `wz`, `Txx`, `Txy`, `Tyy`

This document describes the explicit-Euler time step implemented in
[`src/mms_viscoelastic/numerics/time_stepping.py`](../src/mms_viscoelastic/numerics/time_stepping.py),
exercised by
[`examples/13_case12_euler_step_comparison.py`](../examples/13_case12_euler_step_comparison.py)
and tested in
[`tests/test_case12_euler_step.py`](../tests/test_case12_euler_step.py).

This is **Etapa 9** of the post-FASE-VII replanning
(`docs/planejamento-tecnico.md`, FASE VIII). It ports the right-hand sides of
the `wz`, `Txx`, `Txy` and `Tyy` transport equations and performs one explicit
Euler step starting from the `t=0` MMS Case 12 exact fields, building on the
Poisson solver from Etapa 8
([`docs/case-12-poisson-solver.md`](case-12-poisson-solver.md)).

## 1. Source formulation: `subroutine drv`

The right-hand sides are a verbatim port of `subroutine drv` in
`src/fortran_original/4Order_LPOG_Euler/drvs.f` (the `ty_sim == 1`
formulation), and the explicit update is the Euler branch of
`subroutine solver_euler` in
`src/fortran_original/4Order_LPOG_Euler/euler_rk.f`:

```text
wz  = wz1  + dvz   * dt
Txx = Txx1 + dvTxx * dt
Txy = Txy1 + dvTxy * dt
Tyy = Tyy1 + dvTyy * dt
```

with

```text
dvz = -d(ux*wz)/dx - d(uy*wz)/dy + (beta/Re) * laplacian(wz)
    + d2Txx/dxdy + d2Txy/dy2 - d2Txy/dx2 - d2Tyy/dxdy - tfwz

FTij = 1 + (epsilon*Re*Wi/(1-beta)) * (Txx + Tyy)

dvTxx = -FTij*Txx/Wi - d(ux*Txx)/dx - d(uy*Txx)/dy
      + 2*Txx*dudx + 2*Txy*dudy
      - 2*xi*Txx*dudx - xi*Txy*(dudy + dvdx)
      - (alphaG*Re/(1-beta)) * (Txx^2 + Txy^2)
      + (2*(1-beta)/(Re*Wi)) * dudx + tfTxx/Wi

dvTxy = -FTij*Txy/Wi - d(ux*Txy)/dx - d(uy*Txy)/dy
      + Txx*dvdx + Tyy*dudy
      - xi*Txy*dudx - 0.5*xi*Tyy*(dudy + dvdx)
      - 0.5*xi*Txx*(dudy + dvdx) - xi*Txy*dvdy
      - (alphaG*Re/(1-beta)) * (Txx*Txy + Txy*Tyy)
      + ((1-beta)/(Re*Wi)) * (dudy + dvdx) + tfTxy/Wi

dvTyy = -FTij*Tyy/Wi - d(ux*Tyy)/dx - d(uy*Tyy)/dy
      + 2*Txy*dvdx + 2*Tyy*dvdy
      - xi*Txy*(dvdx + dudy) - 2*xi*Tyy*dvdy
      - (alphaG*Re/(1-beta)) * (Txy^2 + Tyy^2)
      + (2*(1-beta)/(Re*Wi)) * dvdy + tfTyy/Wi
```

`wz_time_derivative` and `stress_time_derivatives`
(`src/mms_viscoelastic/numerics/time_stepping.py`) implement `dvz`, `dvTxx`,
`dvTxy` and `dvTyy` term-by-term from these expressions. `euler_step_wz` and
`euler_step_stress` apply the Euler update.

### Sign convention: `wz = laplacian(psi)`

As documented in
[`docs/case-12-poisson-solver.md`](case-12-poisson-solver.md#2-sign-convention-wz--laplacianpsi),
this package follows `case_12.f`'s convention `wz = laplacian(psi)`. `dvz`
above is the verbatim port of `drv`'s `dvz` under that convention; no sign
flip is introduced.

## 2. Parameter choice

`EulerStepParameters` defaults are:

```text
reynolds = 1.0
wi       = 1.0
beta     = 0.5
epsilon  = 0.0
xi       = 0.0
alpha_g  = 0.0
dt       = 1.0e-5
```

These resolve a discrepancy between two sources in the original Fortran
project:

- `src/fortran_original/4Order_LPOG_Euler/par.nn` sets `Rey=100.0`,
  `Wi=0.05`, `betann=0.5`, `epsylon=1.0`, `xi=0.0`, `alphaG=0.0`,
  `dt=1.0e-05`.
- The thesis figure filenames (e.g.
  `Thesis/figures/Newtonian/..._Betann_1_Re_1_Wi_1_..._Dt_1e-05_..._MMS_13_...`)
  imply `Re=1`, `Wi=1`, `beta=1` (Newtonian limit) for some plotted cases.

`reynolds=1.0`, `wi=1.0`, `beta=0.5`, `epsilon=0.0`, `xi=0.0`, `alpha_g=0.0`
match `Case12Parameters()` and `Case12SourceTermParameters()`, the values
already validated against the Fortran field and source-term references in
FASE V-VII (see
[`docs/case-12-complete-validation-report.md`](case-12-complete-validation-report.md)).
`dt=1.0e-5` matches **both** `par.nn` and the thesis filenames
(`Dt_1e-05`). This is the parameter set used for Etapa 9, and the input CSV
generator (`validation/fortran_reference/case12_euler_step/generate_input_csv.py`)
and the Fortran reference program use the same defaults.

## 3. Discretization

Both the Python implementation and the Fortran reference program
(`validation/fortran_reference/case12_euler_step/case12_euler_step_reference.f90`)
use the **same second-order finite-difference operators** from
[`src/mms_viscoelastic/numerics/derivatives.py`](../src/mms_viscoelastic/numerics/derivatives.py):

- `derivative_x` / `derivative_y`: central differences with one-sided
  second-order edge formulas for the first derivatives.
- `second_derivative_x` / `second_derivative_y`: 3-point central differences
  for the second derivatives, with boundary entries set to zero.
- `laplacian`: the existing 5-point second-order Laplacian, also zero on the
  boundary.
- Mixed derivatives (`d2Txx/dxdy`, `d2Tyy/dxdy`) are computed by composition,
  e.g. `derivative_y(derivative_x(txx, dx), dy)`.

The **4th-order LPOG compact (Padé) operators**
(`derx`/`dery`/`derxx`/`deryy` in `src/fortran_original/4Order_LPOG_Euler/derivs.f`,
solved via the `tridx`/`tridy` Thomas-algorithm routines in `solvers.f`) and
the **Lele compact filter** (`filter_trid_x`/`filter_trid_y` in `loop.f`) are
**not** ported. This is the same deliberate scope reduction documented for the
Poisson solver in Etapa 8
([`docs/case-12-poisson-solver.md`](case-12-poisson-solver.md#simplification-relative-to-the-fortran-reference)).
Porting the 4th-order compact operators and the filter remains deferred.

## 4. Boundary handling and scope

Only **interior grid points** (`[1:-1, 1:-1]` in Python, `i=2..nx-1,
j=2..ny-1` in Fortran) are advanced by the Euler step; boundary values are
copied unchanged from the input fields.

The full Fortran `solver_euler` step also calls, after the Euler update:

```fortran
call source_term(time + dt)
call boundaries_condictions(time + dt)
call loop(2)
```

i.e. it re-evaluates the manufactured source terms at `t = dt`, re-imposes
the Dirichlet boundary conditions at `t = dt`, re-solves the Poisson equation
for `psi` from the updated `wz`, recomputes `ux`/`uy` from `psi`, and applies
the Lele compact filter. **None of this is in scope for Etapa 9** - the
comparison here isolates the RHS-formula port and the Euler arithmetic from
the boundary-condition/Poisson-resolve/filter machinery. Re-imposing boundary
conditions and recomputing `psi`/`ux`/`uy` at `t=dt` is left to the
multi-step time-marching solver (Etapa 10).

## 5. Fortran reference program

`validation/fortran_reference/case12_euler_step/case12_euler_step_reference.f90`
is a standalone program (it does not modify or depend on the original thesis
Fortran files) that:

1. Reads an input CSV with columns `x,y,u,v,wz,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy`,
   generated by `generate_input_csv.py` from the already-validated Python
   `evaluate_case_12_exact` and `evaluate_case_12_source_terms` at `t=0` on a
   `9 x 9` grid. This avoids re-transcribing the large `tfTxx`/`tfTxy`/`tfTyy`
   symbolic expressions into Fortran a second time.
2. Computes `dvz`, `dvTxx`, `dvTxy`, `dvTyy` using the second-order operators
   described above.
3. Advances the interior points by one Euler step of size `dt = 1.0e-5`.
4. Writes an output CSV with columns `x,y,wz,txx,txy,tyy` (post-step values,
   full grid; boundary values equal the input).

See
[`validation/fortran_reference/case12_euler_step/README.md`](../validation/fortran_reference/case12_euler_step/README.md)
for compile and run instructions.

## 6. Observed results

`compare_case12_euler_step_reference` compares each of `wz`, `Txx`, `Txy`,
`Tyy` after one Euler step (`dt = 1.0e-5`, `9 x 9` grid) between Python and
the Fortran reference:

```text
  wz: max_abs_error=4.440892e-16 max_rel_error=2.220446e-16 passed=True
 txx: max_abs_error=4.440892e-16 max_rel_error=1.925336e-16 passed=True
 txy: max_abs_error=4.440892e-16 max_rel_error=2.220443e-16 passed=True
 tyy: max_abs_error=4.440892e-16 max_rel_error=1.925334e-16 passed=True
```

All four components agree to within machine precision (`~4.4e-16`), as
expected since both implementations use identical formulas and discretization
operators applied to identical input fields.
[`tests/test_case12_euler_step.py`](../tests/test_case12_euler_step.py) asserts
`max_abs_error < 1e-12` and `max_rel_error < 1e-12` for all four components.

## 7. Run

```powershell
python validation\fortran_reference\case12_euler_step\generate_input_csv.py `
  validation\reference_data\case12_euler_step_input_9x9.csv 9 9

.\validation\fortran_reference\case12_euler_step\case12_euler_step_reference.exe `
  validation\reference_data\case12_euler_step_input_9x9.csv `
  validation\reference_data\case12_euler_step_output_9x9_dt1e-05.csv `
  9 9 1.0e-5

python examples\13_case12_euler_step_comparison.py
```

## 8. Deferred work

- Etapa 10 (multi-step time marching) re-imposes boundary conditions and
  recomputes `psi`/`ux`/`uy` (via the Poisson solver from Etapa 8) at each
  time step, building on `euler_step_wz` and `euler_step_stress`.
- Porting the 4th-order LPOG compact derivative operators and the Lele
  compact filter remains out of scope, as noted in Section 3.
