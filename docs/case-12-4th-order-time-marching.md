# MMS Case 12 - 4th-order time-marching pipeline (Etapa 17)

This document covers **Etapa 17** of `.claude/PLANEJAMENTO.md`: reapointing
the Case 12 time-marching pipeline from the 2nd-order placeholder (Etapas
8-11) to the real scheme the original Fortran code uses -- the 4th-order
LPOG compact derivative operators (Etapa 13), the Lele compact filter
(Etapa 13), and the FAS multigrid V-cycle Poisson solver -- and rewriting
the four isolated Fortran reference programs accordingly.

## 1. Finding: `band5_poi`/`banbk5_poi` is dead code

Etapa 17 was originally planned (see the Etapa 13 closing note) as "port the
pentadiagonal Poisson solver `band5_poi`/`banbk5_poi` (`solvers.f`) and
reapoint `numerics/poisson.py`". Before implementing that, a repository-wide
search

```text
grep -rn "band5_poi\|banbk5_poi" src/fortran_original/4Order_LPOG_Euler/*.f
```

shows these two subroutines are defined but **never called** anywhere in the
original Fortran source. They are dead code. The Poisson solver actually
invoked by the time-marching loop is `poi_psi4a`
(`call poi_psi4a` inside `subroutine loop`, `main.f:182`): a **FAS (Full
Approximation Scheme) multigrid V-cycle** over a 4th-order compact
(Mehrstellen-type) 9-point discretization, smoothed by line-relaxation SOR
(one tridiagonal solve per grid line, via `tridv`).

Etapa 17 therefore ports `poi_psi4a` and its dependents, not
`band5_poi`/`banbk5_poi`.

## 2. The multigrid Poisson solver

See
[`src/mms_viscoelastic/numerics/poisson_multigrid.py`](../src/mms_viscoelastic/numerics/poisson_multigrid.py)
(`solve_poisson_multigrid_4th_order`) and
[`validation/fortran_reference/case12_poisson_multigrid/README.md`](../validation/fortran_reference/case12_poisson_multigrid/README.md)
for the full algorithm description (the 9-point compact stencil, the FAS
coarse-grid RHS construction, restriction/prolongation). Key points:

- Ports `poi_psi4a`, `verifica_poi`, `function_source`, `sor`, `defect`,
  `restrict`, `intpol`, `correction` and `tridv` verbatim, including the
  exact V-cycle control flow (down-leg smoothing + residual + FAS restrict,
  coarsest-level heavy smoothing, up-leg correction + prolongation +
  smoothing), the `erromax = 1e-6` convergence tolerance and the 81-V-cycle
  cap.
- The line-relaxation sweep in `sor` is **not vectorizable** across grid
  lines: line `i` reads the already-updated line `i-1` and the
  not-yet-updated line `i+1` within the same sweep (Gauss-Seidel zebra
  relaxation, not Jacobi). The Python port loops over lines explicitly to
  reproduce this dependency exactly.
- Validated against an isolated Fortran reference at **exact bit-for-bit
  agreement** (`max_abs_error = 0.0`) on the `33 x 33` grid matching
  `par.nn`'s `imax = jmax = 33`, `msh = 4`, and at 4th-order convergence
  (`Linf` order `~4.0-4.02`, `L2` order `~3.97-4.0`) across `n = 8, 16, 32,
  64` against the Case 12 exact `psi`, both well above the `~1.96-2.01`
  order of the pre-existing 2nd-order direct solver
  (`docs/case-12-poisson-solver.md`).

## 3. Reapointing `time_stepping.py` / `time_stepping_log.py`

`drvs.f`'s `drv` and `drv_Psi` call `derx`/`dery`/`derxx`/`deryy` -- the
4th-order LPOG compact operators from `derivs.f` (already ported in Etapa
13 as `compact_derivative_x/y`/`compact_second_derivative_x/y`), not a
2nd-order scheme. `numerics/time_stepping.py` and
`numerics/time_stepping_log.py` now import these compact operators (aliased
to the module-local names `derivative_x`/`derivative_y`/
`second_derivative_x`/`second_derivative_y`, so the formula code in
`wz_time_derivative`/`stress_time_derivatives`/
`log_conformation_time_derivatives` is unchanged -- only the underlying
discretization is). `laplacian(wz, dx, dy)` is replaced by
`second_derivative_x(wz, dx) + second_derivative_y(wz, dy)` using the
compact operators (no compact "laplacian" helper exists; composition is
equivalent).

This is a behavior change for every function in these two modules -- see
Section 5 for why the existing Fortran references had to be rewritten as a
result.

## 4. Reapointing `time_marching.py`

`subroutine loop(num)` (`main.f:171`) is, in order:

```text
call filter_trid_x
call filter_trid_y
call poi_psi4a
call calc_ux
call calc_uy
```

called right after `boundaries_condictions`, on **every** Euler/RK4 step or
sub-stage (the `num` argument is unused in the body -- a vestige). This
resolves the open question left in Etapa 13 ("should the Lele filter apply
every step, or only at the end?"): every step, between the boundary
re-imposition and the Poisson re-solve.

`case12_time_step` (`numerics/time_marching.py`) now follows this order
exactly:

1. Euler-update the interior of `wz`, `Txx`, `Txy`, `Tyy` (unchanged step,
   now using the compact derivatives per Section 3).
2. Re-impose the exact boundary at `t + dt` (unchanged).
3. Apply the Lele filter, x-pass then y-pass
   (`_apply_lele_filter`/`lele_filter_x`/`lele_filter_y`), to the interior
   of `wz`, `Txx`, `Txy`, `Tyy`. Each Fortran filter pass computes the
   filtered value over the *entire* grid (its one-sided boundary-closure
   stencil needs the boundary values) but only overwrites the interior --
   the boundary stays exactly the value imposed in step 2. The y-pass reads
   the already x-filtered field, so the two passes are sequential, not
   commutative or mergeable.
4. Re-solve `laplacian(psi) = wz` with
   `solve_poisson_multigrid_4th_order` (Section 2), using the filtered `wz`
   and the exact `psi` boundary/initial guess.
5. Recompute `u = d(psi)/dy`, `v = -d(psi)/dx` with the compact derivatives
   (`compact_derivative_y`/`compact_derivative_x`), matching `calc_ux`/
   `calc_uy` (`loop.f`).

The `wz`-boundary deviation from `boundarie_vorticity_cavity_flow.f`
(documented since Etapa 10) is unchanged by Etapa 17.

## 5. Rewriting the Fortran reference programs

Because Sections 3-4 change the actual numerical behavior of
`time_stepping.py`/`time_stepping_log.py`/`time_marching.py`, the five
existing isolated Fortran references that validate them against the
2nd-order placeholder no longer matched -- they had to be rewritten to use
the same real scheme, not just have their tolerances loosened:

| Reference | Change |
|---|---|
| `case12_euler_step` | `derivative_x`/`derivative_y`/`second_derivative_x`/`second_derivative_y`/`laplacian` (same names/signatures) reimplemented via the compact operators (`coef_4a`/`coeff_4a`/`rhsx_4a`/`rhsxx_4a`/`rhsy_4a`/`rhsyy_4a`/`tridx`/`tridy`, copied from `case_compact_derivatives`); no other call sites change. |
| `case12_log_euler_step` | Same derivative-subroutine swap as above (no filter/Poisson involved at this stage). |
| `case12_time_step` | Same derivative swap, **plus** the dense 2nd-order Poisson assembly + Gaussian-elimination solve replaced by `apply_lele_filter` (x then y) on `wz`/`Txx`/`Txy`/`Tyy`, followed by `poi_psi4a` (copied verbatim from `case12_poisson_multigrid`). |
| `case12_multistep` | Same as `case12_time_step`, applied inside the existing `do step = 1, n_steps` loop. |
| `case12_full_march` | Same as `case12_time_step`/`case12_multistep`, applied inside the existing `do t_step = 1, max_steps` convergence loop. See Section 6 for what changed in its long-horizon behavior. |

All five were recompiled and their reference CSVs regenerated. Every
existing test continues to pass at the same machine-precision tolerances
(`< 1e-12`) as before -- the comparison is now between two 4th-order
implementations instead of two 2nd-order ones:

```text
case12_euler_step:      wz/Txx/Txy/Tyy        max_abs ~4.4e-16
case12_log_euler_step:  wz/Txx/Txy/Tyy        max_abs ~4.4e-16 to ~4.0e-15
case12_time_step:       wz/Txx/Txy/Tyy/psi/u/v max_abs ~2.8e-17 to ~4.4e-16
case12_multistep (n=10): same components       max_abs ~4.2e-17 to ~9.4e-16
case12_poisson_multigrid: psi                  max_abs = 0.0 (exact)
```

## 6. `case12_full_march`: the filter does not prevent divergence either

`docs/case-12-full-march-reference.md` (Etapa 12) found that the 2nd-order,
no-filter discretization **diverges** on a `9 x 9` grid at step ~455064
(`time ~= 4.5506`), and hypothesized that the real, dissipative
discretization (4th-order + Lele filter + multigrid) was "not expected to
share this instability."

Running the rewritten `case12_full_march` with the real scheme on the same
`9 x 9` grid (`dt = 1e-5`) **refutes that hypothesis**: it also diverges,
in fact *earlier* -- the multigrid solver itself fails to converge within
its 81-V-cycle cap (`error stop 'Vcycles = 82'`) by around step ~238000
(`time ~= 2.38`), with the field maxima already growing explosively
(`max|Txx|` from `O(10)` to `O(10^3)` over the preceding ~16000 steps,
the same qualitative blow-up pattern as the 2nd-order case, just on a
different timescale). The Lele filter removes grid-scale (odd-even) noise;
it does not address the underlying explicit-Euler-at-fixed-`dt`-on-a-coarse-
grid instability.

On a `17 x 17` grid, the real scheme runs the full `500000` steps
(`time = 5.0`) without diverging (`diverged = F`), consistent with the
original finding that finer grids are more stable -- but `emv` does not
reach `solver_euler`'s `1e-7` convergence threshold either (`final emv =
9.4e-05`), so "converges to discrete steady state" is not reached within a
practical horizon at this resolution.

**Conclusion**: the instability documented in
`docs/case-12-full-march-reference.md` is a property of explicit Euler at
`dt = 1e-5` on these grid resolutions, essentially independent of whether
the spatial discretization is the 2nd-order placeholder or the real
4th-order-compact + Lele-filter + multigrid-Poisson scheme. Etapa 13's
re-scoping decision (validate `case12_multistep_march` over a stable,
bounded horizon, not to convergence or divergence) remains the right call
for both schemes. No automated test depends on `case12_full_march`'s
long-horizon CSVs (only this documentation references them), so they were
not regenerated at full length; the exploratory runs above were not saved
as committed artifacts.

## 7. Run

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_poisson_multigrid\case12_poisson_multigrid_reference.f90 `
  -o validation\fortran_reference\case12_poisson_multigrid\case12_poisson_multigrid_reference.exe
.\validation\fortran_reference\case12_poisson_multigrid\case12_poisson_multigrid_reference.exe `
  validation\reference_data\case12_poisson_multigrid_input_33x33.csv `
  validation\reference_data\case12_poisson_multigrid_output_33x33.csv 33 33 4

gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_euler_step\case12_euler_step_reference.f90 `
  -o validation\fortran_reference\case12_euler_step\case12_euler_step_reference.exe
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_time_step\case12_time_step_reference.f90 `
  -o validation\fortran_reference\case12_time_step\case12_time_step_reference.exe
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_multistep\case12_multistep_reference.f90 `
  -o validation\fortran_reference\case12_multistep\case12_multistep_reference.exe

python -m pytest tests\test_poisson_multigrid.py `
  tests\test_case12_poisson_multigrid_fortran_reference.py `
  tests\test_case12_euler_step.py `
  tests\test_case12_log_conformation_fortran_reference.py `
  tests\test_case12_time_step.py `
  tests\test_case12_multistep_march.py
```

## 8. Deferred work

- `case12_full_march`'s long-horizon trajectories (Section 6) were explored
  but not committed as regenerated CSVs/history files, since no test reads
  them.
- The RK4 variants (`solver_rk4`, `solver_log_rk4`) still use the 2nd-order
  and log-conformation-Euler-only ports respectively, unaffected by this
  Etapa.
- `cc_Psi` (log-conformation boundary re-projection) remains unported
  (Etapa 15 scope).
