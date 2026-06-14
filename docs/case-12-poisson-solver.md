# MMS Case 12 - Poisson solver for `psi`

This document describes the direct Poisson solver implemented in
[`src/mms_viscoelastic/numerics/poisson.py`](../src/mms_viscoelastic/numerics/poisson.py),
exercised by
[`examples/12_case12_poisson_solver.py`](../examples/12_case12_poisson_solver.py)
and tested in
[`tests/test_case12_poisson_solver.py`](../tests/test_case12_poisson_solver.py).

This is **Etapa 8** of the post-FASE-VII replanning
(`docs/planejamento-tecnico.md`, FASE VIII). It provides the missing building
block for a Python time-marching solver: recomputing the stream function
`psi` from the vorticity `wz` via `laplacian(psi) = wz`.

## 1. Why this is needed

[`docs/case-12-solver-comparison.md`](case-12-solver-comparison.md) (FASE VII)
compares the Python `evaluate_case_12_exact` against the Fortran reference
only at `t = 0`, because the Python package had no discretized solver. A
time-marching loop needs, at every step, a way to recover `psi` from the
updated vorticity field `wz` - this is exactly the 2D Poisson problem with
Dirichlet boundary conditions:

```text
laplacian(psi) = wz   in the interior
psi = psi_boundary     on the domain boundary
```

## 2. Sign convention: `wz = laplacian(psi)`

[`docs/case-12-complete-validation-report.md`](case-12-complete-validation-report.md#12-a-note-on-the-sign-of-wz)
documents that, in this package (following `case_12.f`), `wz` is defined as

```text
wz = d(u)/dy - d(v)/dx = laplacian(psi)
```

which is the **negative** of the vorticity `omega_z` in the README's general
formulation. The Poisson solver in this module follows the same convention:
it solves `laplacian(phi) = source` with `source = wz`, **not** `-wz`. Any
future code that calls `solve_poisson_dirichlet` to recompute `psi` from `wz`
must pass `wz` directly as `source`.

## 3. Discretization

`solve_poisson_dirichlet(source, grid, boundary_values)` discretizes
`laplacian(phi) = source` with the standard second-order, five-point stencil

```text
(phi[i, j-1] - 2*phi[i, j] + phi[i, j+1]) / dx^2
  + (phi[i-1, j] - 2*phi[i, j] + phi[i+1, j]) / dy^2 = source[i, j]
```

for every interior point `(i, j)`. Known boundary values from
`boundary_values` are moved to the right-hand side for boundary-adjacent
interior points. The resulting sparse linear system (Kronecker sum of two 1D
second-difference matrices, built with `scipy.sparse`) is solved directly
with `scipy.sparse.linalg.spsolve`. The function returns a full
`(grid.ny, grid.nx)` array: the boundary entries are copied verbatim from
`boundary_values`, and the interior entries are the direct-solve result.

### Simplification relative to the Fortran reference

The original Fortran solver
(`src/fortran_original/4Order_LPOG_Euler/poisson.f`) uses a multigrid V-cycle
with SOR smoothing (`poi_psi4a`, 4 grid levels, convergence tolerance
`erromax = 1.d-06`). For the manufactured-solution verification grids used
here (up to ~65 x 65 points), a sparse direct solve gives the same
second-order discretization and the same solution (to round-off), without the
complexity of porting the multigrid algorithm. The multigrid solver itself is
**not** ported; this is documented here as a deliberate scope reduction (see
`.claude/PLANEJAMENTO.md`, Etapa 8).

## 4. MMS verification

`run_case12_poisson_convergence_study(grid_sizes, t=0.0, parameters=None)`
verifies the solver against case 12:

1. For each grid `n` (cells per direction, `nx = ny = n + 1` points), evaluate
   the exact fields `wz_exact` and `psi_exact` via `evaluate_case_12_exact`.
2. Solve `laplacian(psi_numerical) = wz_exact` with `psi_exact` as the
   Dirichlet boundary condition.
3. Compute the `L2` and `Linf` errors between `psi_numerical` and
   `psi_exact` over the whole grid, and estimate the order of accuracy via a
   log-log fit across grids.

### Observed results

```text
    n            h         L2_psi       Linf_psi
    8 1.250000e-01   1.900267e-03   4.981034e-03
   16 6.250000e-02   4.962755e-04   1.229010e-03
   32 3.125000e-02   1.273984e-04   3.062456e-04
   64 1.562500e-02   3.231234e-05   7.649858e-05

Estimated order of accuracy (log-log fit across all grids):
  psi: L2 order = 1.9596, Linf order = 2.0079
```

Both `L2` and `Linf` errors decrease monotonically and the estimated order is
close to `2`, consistent with the second-order discretization. This matches
the order observed for the stream-function identities in
[`docs/case-12-convergence-verification.md`](case-12-convergence-verification.md).

## 5. Acceptance criteria

[`tests/test_case12_poisson_solver.py`](../tests/test_case12_poisson_solver.py)
checks, for grids `n = 8, 16, 32, 64`:

- `solve_poisson_dirichlet` exactly reproduces a known quadratic
  (`phi = x^2 + y^2`, `laplacian(phi) = 4`) to within `1e-10`.
- shape-mismatch inputs raise `ValueError`.
- `L2` errors of `psi_numerical` vs `psi_exact` decrease monotonically with
  refinement.
- `L2_order >= 1.8` and `Linf_order >= 1.8`.
- at `n = 32` (33 x 33 points), `max(|psi_numerical - psi_exact|) < 1e-3`.

## 6. Run

```powershell
python examples\12_case12_poisson_solver.py
python examples\12_case12_poisson_solver.py --grid-sizes 8 16 32 64 128
```

## 7. Deferred work

- Etapa 9 (single Euler time step) and Etapa 10 (one full time-marching step)
  build on this solver to compute `psi` from the evolving `wz` field. See
  [`docs/case-12-euler-step-comparison.md`](case-12-euler-step-comparison.md)
  and
  [`docs/case-12-time-marching-comparison.md`](case-12-time-marching-comparison.md).
