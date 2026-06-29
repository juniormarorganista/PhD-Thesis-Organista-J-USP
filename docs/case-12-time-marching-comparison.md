# MMS Case 12 - One full time-marching step

This document describes the time-marching step implemented in
[`src/mms_viscoelastic/numerics/time_marching.py`](../src/mms_viscoelastic/numerics/time_marching.py),
exercised by
[`examples/14_case12_time_step_comparison.py`](../examples/14_case12_time_step_comparison.py)
and tested in
[`tests/test_case12_time_step.py`](../tests/test_case12_time_step.py).

This is **Etapa 10** of the post-FASE-VII replanning
(`docs/planejamento-tecnico.md`, FASE VIII). It combines the explicit Euler
step from Etapa 9
([`docs/case-12-euler-step-comparison.md`](case-12-euler-step-comparison.md))
with the Poisson solver from Etapa 8
([`docs/case-12-poisson-solver.md`](case-12-poisson-solver.md)) into one full
time step: Euler update, boundary-condition re-imposition at `t+dt`, Poisson
re-solve for `psi`, and recomputation of `u`/`v` from the updated `psi`.

> **Etapa 17 update.** Section 3 below described a 2nd-order discretization
> and a direct sparse Poisson solve as a deliberate scope reduction. As of
> Etapa 17, `case12_time_step` uses the real pipeline: 4th-order LPOG
> compact derivatives, the Lele compact filter applied between the
> boundary re-imposition and the Poisson re-solve, and the FAS multigrid
> V-cycle Poisson solver. The Fortran reference program and its CSV were
> rewritten/regenerated to match. See
> [`docs/case-12-4th-order-time-marching.md`](case-12-4th-order-time-marching.md)
> for the full story. Sections 1-2 and 4-7 below are otherwise unchanged.

## 1. Sequence ported

`case12_time_step` is a verbatim-in-spirit port of the sequence that follows
the Euler update in `subroutine solver_euler`
(`src/fortran_original/4Order_LPOG_Euler/euler_rk.f`):

```fortran
call source_term(time + dt)
call boundaries_condictions(time + dt)
call loop(2)
```

Concretely, given the field state `(u, v, wz, psi, Txx, Txy, Tyy)` at time
`t`:

1. Evaluate the manufactured source terms `tfwz`, `tfTxx`, `tfTxy`, `tfTyy` at
   `t` (`evaluate_case_12_source_terms`).
2. Advance `wz`, `Txx`, `Txy`, `Tyy` at the interior points by one explicit
   Euler step (`euler_step_wz`, `euler_step_stress`, Etapa 9).
3. Re-impose Dirichlet boundary conditions at `t + dt` for `wz`, `Txx`, `Txy`,
   `Tyy`, `psi`, `u` and `v` from `evaluate_case_12_exact(grid, t + dt)`.
4. Re-solve `laplacian(psi) = wz` for the interior of `psi`
   (`solve_poisson_dirichlet`, Etapa 8), using the `psi` boundary values from
   step 3.
5. Recompute `u = d(psi)/dy` and `v = -d(psi)/dx` at the interior points from
   the updated `psi` (the same identities verified in
   [`docs/case-12-convergence-verification.md`](case-12-convergence-verification.md)),
   then re-impose the `u`/`v` boundary values from step 3.

`case12_time_step` returns a new `TimeMarchingState` at `t + dt`. Repeated
calls (passing `t = t + dt` each time) perform multi-step time marching;
chaining steps is left to future work (see Section 6).

## 2. Boundary conditions: `boundaries_condictions_case_12`

`subroutine boundaries_condictions_case_12(t)`
(`src/fortran_original/4Order_LPOG_Euler/case_12.f`, lines 360-468) re-imposes
the exact manufactured solution on the domain boundary for `u`, `v`, `psi`,
`Txx`, `Txy` and `Tyy` at every time step, using the same formulas as
`exact_solutions_case_12(t)` (ported in `evaluate_case_12_exact`). This part
is a direct port: step 3 above simply overwrites the outer boundary rows and
columns of `wz`, `Txx`, `Txy`, `Tyy` and `psi` (and, after the Poisson
re-solve, of `u` and `v`) with `evaluate_case_12_exact(grid, t + dt)`.

### `wz` boundary condition: a deliberate simplification

The Fortran subroutine does **not** set `wz` from the exact solution.
Instead, it calls `boundarie_vorticity_cavity_flow`
(`src/fortran_original/4Order_LPOG_Euler/boundarie_vorticity_cavity_flow.f`),
which derives the `wz` boundary values from one-sided finite-difference
stencils on `psi` (bottom/left/right boundaries) and `u` (top boundary), e.g.

```fortran
dy2 = 18.d0*dy*dy
wz(i,1) = (-85*psi(i,1) + 108*psi(i,2) - 27*psi(i,3) + 4*psi(i,4)) / dy2
```

These stencils are 4th-order one-sided formulas tied to the LPOG compact
discretization (`derx`/`dery`/`derxx`/`deryy` in `derivs.f`), which Etapas 8
and 9 already left out of scope.

`case12_time_step` instead re-imposes the `wz` boundary directly from
`evaluate_case_12_exact(grid, t + dt).wz`, the same way as the other boundary
fields. Both choices converge to the same exact boundary value as the grid is
refined; the Fortran one-sided stencil is simply a discrete approximation of
it consistent with the 4th-order interior operators. Using the exact value is
consistent with the 2nd-order scope reduction already documented for the
Poisson solver (Etapa 8) and the Euler step (Etapa 9). **This is a deliberate
deviation from the Fortran reference**, and the Fortran reference program used
for comparison (Section 4) makes the same choice, so the comparison in
Section 5 isolates the rest of the time-step formulation.

## 3. Discretization

As in Etapas 8 and 9, both the Python implementation and the Fortran
reference program use the **same second-order finite-difference operators**
from
[`src/mms_viscoelastic/numerics/derivatives.py`](../src/mms_viscoelastic/numerics/derivatives.py)
for the Euler step right-hand sides and for recomputing `u`/`v` from `psi`.

The Poisson re-solve in Python uses
`solve_poisson_dirichlet` (`scipy.sparse` + `spsolve`, Etapa 8). The Fortran
reference assembles the same 5-point second-order Laplacian system (identical
ordering and boundary-term construction) and solves it with a dense Gaussian
elimination with partial pivoting -- both are exact solves of the same linear
system, so they agree to round-off.

The 4th-order LPOG compact operators and the Lele compact filter remain out of
scope, as in Etapas 8 and 9.

## 4. Fortran reference program

[`validation/fortran_reference/case12_time_step/case12_time_step_reference.f90`](../validation/fortran_reference/case12_time_step/case12_time_step_reference.f90)
is a standalone program (it does not modify or depend on the original thesis
Fortran files) that:

1. Reads an input CSV with columns
   `x,y,u,v,wz,psi,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy`, generated by
   [`generate_input_csv.py`](../validation/fortran_reference/case12_time_step/generate_input_csv.py)
   from the already-validated Python `evaluate_case_12_exact` and
   `evaluate_case_12_source_terms` at `t=0` on a `9 x 9` grid (extending the
   Etapa 9 input CSV with a `psi` column).
2. Performs the explicit Euler step (identical to Etapa 9).
3. Evaluates the exact manufactured fields at `t = dt` (the same formulas
   ported for `validation/fortran_reference/case12_fields/case12_fields_reference.f90`)
   and uses them to re-impose the boundary conditions on `wz`, `Txx`, `Txy`,
   `Tyy` and `psi`.
4. Assembles and solves the interior Poisson system for `psi` by dense
   Gaussian elimination.
5. Recomputes `u`/`v` from the updated `psi` and re-imposes their boundary
   values.
6. Writes an output CSV with columns `x,y,wz,txx,txy,tyy,psi,u,v` (post-step
   values, full grid).

See
[`validation/fortran_reference/case12_time_step/README.md`](../validation/fortran_reference/case12_time_step/README.md)
for compile and run instructions.

## 5. Observed results

`compare_case12_time_step_reference` compares each of `wz`, `Txx`, `Txy`,
`Tyy`, `psi`, `u`, `v` after one full time-marching step (`dt = 1.0e-5`,
`9 x 9` grid) between Python and the Fortran reference:

```text
  wz: max_abs_error=4.440892e-16 max_rel_error=2.220446e-16 passed=True
 txx: max_abs_error=4.440892e-16 max_rel_error=1.925336e-16 passed=True
 txy: max_abs_error=4.440892e-16 max_rel_error=2.220443e-16 passed=True
 tyy: max_abs_error=4.440892e-16 max_rel_error=1.925334e-16 passed=True
 psi: max_abs_error=6.938894e-17 max_rel_error=6.938894e-17 passed=True
   u: max_abs_error=2.220446e-16 max_rel_error=2.220446e-16 passed=True
   v: max_abs_error=2.775558e-16 max_rel_error=2.775558e-16 passed=True
```

All seven components agree to within machine precision (`~4.4e-16`), as
expected since both implementations use identical formulas, the same
second-order operators, and exact solves of the same discrete Poisson system.
[`tests/test_case12_time_step.py`](../tests/test_case12_time_step.py) asserts
`max_abs_error < 1e-12` and `max_rel_error < 1e-12` for all seven components.

## 6. Run

```powershell
python validation\fortran_reference\case12_time_step\generate_input_csv.py `
  validation\reference_data\case12_time_step_input_9x9.csv 9 9

.\validation\fortran_reference\case12_time_step\case12_time_step_reference.exe `
  validation\reference_data\case12_time_step_input_9x9.csv `
  validation\reference_data\case12_time_step_output_9x9_dt1e-05.csv `
  9 9 1.0e-5

python examples\14_case12_time_step_comparison.py
```

## 7. Deferred work

- **Multi-step time marching**: implemented in Etapa 11
  ([`docs/case-12-multistep-marching-comparison.md`](case-12-multistep-marching-comparison.md)).
  `case12_multistep_march` chains `case12_time_step` over `n_steps`
  increments of `dt`, and is compared against a multi-step Fortran reference.
- Porting the 4th-order LPOG compact derivative operators, the
  `boundarie_vorticity_cavity_flow` one-sided `wz`-boundary stencil, and the
  Lele compact filter remains out of scope, as noted in Sections 2 and 3.
