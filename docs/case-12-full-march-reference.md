# MMS Case 12 - Full time march reference (Etapa 12)

This document is **Etapa 12** of `docs/planejamento-tecnico.md` (FASE IX). It
describes
[`validation/fortran_reference/case12_full_march/case12_full_march_reference.f90`](../validation/fortran_reference/case12_full_march/case12_full_march_reference.f90),
a standalone Fortran program that extends the Etapa 10/11 time-marching step
(explicit Euler + boundary re-imposition + Poisson re-solve + `u`/`v`
recompute, second-order operators) into a loop that mirrors the *stopping
criterion* of `solver_euler`
(`src/fortran_original/4Order_LPOG_Euler/euler_rk.f`, lines 326-417):

```text
emv = max(errowz, erroTxx, erroTxy, erroTyy)
if (emv < 1.0e-7 .and. time > 0.01) exit  ! converged to discrete steady state
```

with `errowz = max|wz_new - wz_old| / max(max|wz_new|, 1)` (and analogously
for `Txx`, `Txy`, `Tyy`).

## 1. Summary of the finding

**The simplified second-order reference does not satisfy `solver_euler`'s
convergence criterion. On a `9 x 9` grid it diverges (overflows) at step
455064 (`time ~= 4.5506`) instead of converging.**

This is a genuine numerical instability of the *explicit Euler + 2nd-order
central differences, without the Lele compact filter*, at `dt = 1e-5` -- not
a translation bug. Etapas 9-11 already validated that this reference and
`case12_multistep_march` agree to `~1e-16` for short runs (10 steps), so the
per-step formulas are correct; the instability is a property of repeating
those (correct) formulas hundreds of thousands of times without the
dissipation that the real solver's Lele filter (`loop(2)`'s
`filter_trid_x`/`filter_trid_y`, `alphaf = 0.48`) provides.

This finding re-scopes Etapas 13-15 (Section 5).

## 2. `9 x 9` trajectory

Run: `case12_full_march_reference.exe ... 9 9 1.0e-5 history.csv 460000 50.0 1 200`.
Full history in
[`validation/reference_data/case12_full_march_history_9x9.csv`](../validation/reference_data/case12_full_march_history_9x9.csv)
(2276 rows, every 200th step), fine-grained tail in
[`validation/reference_data/case12_full_march_history_9x9_tail.csv`](../validation/reference_data/case12_full_march_history_9x9_tail.csv)
(steps 454900-455064, every step). Final state in
[`validation/reference_data/case12_full_march_output_9x9.csv`](../validation/reference_data/case12_full_march_output_9x9.csv).

| step | time | emv | max\|wz\| | max\|Txx\| | max\|psi\| |
|---|---|---|---|---|---|
| 1 | 1.0e-5 | 5.99e-06 | 4.788 | 2.914 | 0.1493 |
| 30001 | 0.30001 | 2.63e-06 | 5.112 | 3.018 | 0.1603 |
| 100001 | 1.00001 | 3.87e-06 | 5.218 | 3.049 | 0.1619 |
| 200001 | 2.00001 | 1.24e-05 | 5.675 | 3.100 | 0.1683 |
| 300001 | 3.00001 | 3.05e-05 | 5.386 | 9.649 | 0.1806 |
| 400001 | 4.00001 | 2.73e-05 | 10.33 | 15.26 | 0.2140 |
| 440001 | 4.40001 | 1.08e-04 | 23.11 | 17.58 | 0.3664 |
| 450001 | 4.50001 | 4.47e-04 | 76.10 | 80.97 | 0.6531 |
| 454401 | 4.54401 | 2.67e-03 | 938.2 | 614.5 | 5.080 |
| 455001 | 4.55001 | 4.15e-02 | 15629 | 81318 | 79.25 |
| 455064 | 4.55064 | NaN -> `diverged=T` | Infinity | 6.0e+230 | 2.3e+221 |

Observations:

- `emv` starts at `5.99e-6`, *decreases* to a minimum of `~2.6e-6` around
  `time ~= 0.3-0.4`, then grows (with minor non-monotonicity) for the rest of
  the run. `solver_euler`'s `1e-7` threshold is **never** reached -- the
  minimum `emv` over the whole trajectory (`~2.6e-6`) is already more than an
  order of magnitude above it.
- `max|wz|` starts at `4.788`, matching the manufactured solution's amplitude
  (`wz_case12_default` is `O(1)` on `[0,1]x[0,1]`). It grows slowly to `~76`
  by `time = 4.5` (a 16x amplification), then **explosively** by factors of
  `~10x` per `0.001` time units over the last `0.0006` time units, until
  `wz`, `Txx`, `Txy`, `Tyy` and `psi` all overflow to `+-Infinity`/`NaN` at
  step 455064.
- The growth is smooth and monotonic in the tail (see the tail CSV): this is
  consistent with a single unstable eigenmode of the discrete one-step
  operator with `|1 + dt*lambda| > 1` slowly dominating the solution, i.e. a
  textbook linear numerical instability, not a one-off glitch.

## 3. `17 x 17` trajectory (stable horizon)

Run: `case12_full_march_reference.exe ... 17 17 1.0e-5 history.csv 150000 50.0 1 500`
(`time` up to `1.5`, well short of any divergence). Full history in
[`validation/reference_data/case12_full_march_history_17x17.csv`](../validation/reference_data/case12_full_march_history_17x17.csv).

| step | time | emv | max\|wz\| | max\|Txx\| |
|---|---|---|---|---|
| 1 | 1.0e-5 | 1.54e-06 | 4.837 | 2.914 |
| 50001 | 0.50001 | 1.21e-06 | 4.925 | 2.939 |
| 100001 | 1.00001 | 1.61e-06 | 4.912 | 2.941 |
| 149501 | 1.49501 | 2.52e-06 | 4.903 | 2.942 |

`converged = F`, `diverged = F`, `final emv = 2.54e-06` at `time = 1.5`.

Over the same time interval, `17 x 17`'s `emv` (`~1.2e-6` to `~2.5e-6`) is
roughly half of `9 x 9`'s (`~2.6e-6` to `~6.8e-6`, interpolating Section 2),
and `max|wz|` stays within `2%` of its initial value (vs. `9 x 9`'s `~9%`
drift by `time = 1`). This is consistent with the instability being driven by
the discretization error of the simplified scheme (which decreases with grid
refinement): finer grids grow more slowly, but nothing in this run rules out
eventual divergence on `17 x 17` too at a later time.

## 4. Performance

| grid | steps run | wall time | s/step | notes |
|---|---|---|---|---|
| `9 x 9` (`n=49` interior unknowns) | 455064 (to divergence) | 11.98 s | 2.64e-5 | dense Gaussian elimination, `O(n^3)` per step |
| `17 x 17` (`n=225` interior unknowns) | 150000 (`time=1.5`, stable) | 257.0 s | 1.71e-3 | same solver |

The observed `s/step` ratio (`~65x`) is in the same order of magnitude as the
`O(n^3)` prediction `(225/49)^3 ~= 97x`, with the difference attributable to
lower-order terms of the dense solve at these small sizes.

Extrapolating `17 x 17` to `9 x 9`'s divergence horizon (`455064` steps) would
cost `~3200 s` (`~53 min`) in this Fortran reference alone. A Python port using
`scipy.sparse.linalg.spsolve` with a full matrix rebuild every step (as
`solve_poisson_dirichlet` currently does) would be substantially slower per
step. **Etapa 13 should therefore target a bounded number of steps (Section
5), not a march to divergence or to `solver_euler`'s convergence criterion.**

## 5. Re-scoping Etapas 13-15

The FASE IX plan (`docs/planejamento-tecnico.md`) framed the
"self-consistency axis" (Etapa 13) as validating `case12_multistep_march`
against this reference *at the point where `solver_euler`'s convergence
criterion is met*. Section 1 shows that point does not exist for this
simplified discretization -- the reference diverges instead.

**Revised scope for Etapa 13**: validate `case12_multistep_march` against
this reference over a **stable, bounded horizon** chosen well before the
instability dominates, e.g. `n_steps` up to `100000-200000` (`time` up to
`1.0-2.0`), on `9 x 9` and `17 x 17`. This:

- still delivers "larger grids" and "much longer horizon than Etapa 11's
  10 steps" (a `10000x`-`20000x` increase), which was the actual gap left by
  FASE VIII;
- stays in the region where both `max|wz|` and `emv` are close to their
  initial values (Sections 2-3), so a rounding-error-level comparison
  (`~1e-15`, as in Etapas 9-11) is meaningful;
- avoids the open question of whether Python (`numpy`/`scipy`) and Fortran
  (`gfortran`) accumulate rounding differences identically over hundreds of
  thousands of steps of a system with at least one slowly-growing unstable
  mode -- a march to the `9 x 9` divergence point (`455064` steps) would risk
  amplifying `~1e-16` rounding differences between the two implementations
  into `O(1)` disagreements well before either one "diverges", which would be
  meaningless as a validation result.

**Etapas 14-15 (cross-validation axis) are unaffected**: they compare against
the *real, unmodified* `solver_euler` (4th-order LPOG operators + Lele
compact filter + multigrid Poisson), which is a different, dissipative
discretization not expected to share this instability. `solver_euler`'s own
`emv < 1e-7` convergence criterion remains the right stopping condition for
Etapas 14-15.

## 6. Reproduce

```powershell
python validation\fortran_reference\case12_time_step\generate_input_csv.py `
  validation\reference_data\case12_time_step_input_9x9.csv 9 9
python validation\fortran_reference\case12_time_step\generate_input_csv.py `
  validation\reference_data\case12_time_step_input_17x17.csv 17 17

gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_full_march\case12_full_march_reference.f90 `
  -o validation\fortran_reference\case12_full_march\case12_full_march_reference.exe

.\validation\fortran_reference\case12_full_march\case12_full_march_reference.exe `
  validation\reference_data\case12_time_step_input_9x9.csv `
  validation\reference_data\case12_full_march_output_9x9.csv `
  9 9 1.0e-5 validation\reference_data\case12_full_march_history_9x9.csv 460000 50.0 1 200

.\validation\fortran_reference\case12_full_march\case12_full_march_reference.exe `
  validation\reference_data\case12_time_step_input_17x17.csv `
  validation\reference_data\case12_full_march_output_17x17.csv `
  17 17 1.0e-5 validation\reference_data\case12_full_march_history_17x17.csv 150000 50.0 1 500
```

See
[`validation/fortran_reference/case12_full_march/README.md`](../validation/fortran_reference/case12_full_march/README.md)
for the program's CLI and CSV schemas.
