# MMS Case 12 - grid-convergence verification

This document describes the verification implemented in
[`src/mms_viscoelastic/numerics/convergence.py`](../src/mms_viscoelastic/numerics/convergence.py),
exercised by
[`examples/09_case12_convergence_verification.py`](../examples/09_case12_convergence_verification.py)
and tested in
[`tests/test_case12_convergence.py`](../tests/test_case12_convergence.py).

The goal is to confirm, by the Method of Manufactured Solutions (MMS), that
the finite-difference operators in
[`numerics/derivatives.py`](../src/mms_viscoelastic/numerics/derivatives.py)
reproduce known exact differential relations of case 12 at the expected
order of accuracy `p = 2`.

## 1. Exact differential identities satisfied by case 12

The exact fields of case 12,
[`evaluate_case_12_exact`](../src/mms_viscoelastic/manufactured_solutions/case_12.py),
are built from a stream function `psi(x, y, t)`:

```text
psi = x(x - 1) * y * sin(pi x) * sin(pi y) / exp(at * t)
```

together with `u`, `v` and `wz` given by closed-form trigonometric/polynomial
expressions (see the source file for the full formulas). These fields are
not independent: they satisfy the classical 2D stream-function identities

```text
u  = d(psi)/dy                      (1)
v  = -d(psi)/dx                     (2)
wz = d(u)/dy - d(v)/dx              (3)
```

Identity (3) follows from (1)-(2) since

```text
d(u)/dy - d(v)/dx = d^2(psi)/dy^2 + d^2(psi)/dx^2 = laplacian(psi)
```

so `wz = laplacian(psi)` is an equivalent statement of (3).

**Verification.** Identities (1)-(3) were checked **symbolically** with
SymPy by differentiating the closed-form expressions for `u`, `v`, `psi` and
`wz` from `evaluate_case_12_exact` and simplifying

```text
u  - d(psi)/dy            -> 0
v  + d(psi)/dx             -> 0
wz - (d(u)/dy - d(v)/dx)   -> 0
```

All three expressions simplify to exactly `0` for arbitrary `x`, `y`, `t`,
`at`, confirming that (1)-(3) are exact algebraic identities of the
manufactured solution, not approximations. This makes them a clean target
for verifying the *discrete* operators: any deviation between the
finite-difference estimate of (1)-(3) and the exact field is purely
discretization error, with no modeling error contaminating the result.

> Note: identity (3) was checked in the form `wz = d(u)/dy - d(v)/dx`
> (equivalently `wz = laplacian(psi)`), **not** `wz = -laplacian(psi)`. The
> sign matters: case 12 follows the convention in which the vorticity field
> equals `+laplacian(psi)`, opposite to the more common
> `omega = curl(u) = dv/dx - du/dy` convention. Using the wrong sign produces
> an `O(1)` error that does not decrease under refinement.

## 2. Discrete operators and their truncation error

`derivative_x` and `derivative_y`
([`numerics/derivatives.py`](../src/mms_viscoelastic/numerics/derivatives.py))
approximate `d f/dx` (respectively `d f/dy`) using:

- a **central difference** `(f[i+1] - f[i-1]) / (2h)` at interior points;
- a **one-sided three-point** stencil `(-3 f[0] + 4 f[1] - f[2]) / (2h)` (and
  its mirror) at the domain boundaries.

By Taylor expansion, both stencils have leading-order truncation error
`O(h^2)`:

```text
(f(x+h) - f(x-h)) / (2h)              = f'(x) + (h^2/6) f'''(x) + O(h^4)
(-3 f(x) + 4 f(x+h) - f(x+2h)) / (2h) = f'(x) - (h^2/3) f'''(x) + O(h^3)
```

Both leading error terms are `O(h^2)`, with the one-sided stencil carrying a
larger constant (`1/3` vs. `1/6`).

(see e.g. LeVeque, *Finite Difference Methods for Ordinary and Partial
Differential Equations*, SIAM 2007, Ch. 1, for both stencils). Therefore the
discrete relations

```text
u_num  = derivative_y(psi)
v_num  = -derivative_x(psi)
wz_num = derivative_y(u) - derivative_x(v)
```

should converge to the exact fields `u`, `v`, `wz` with error `O(h^2)`,
i.e. order `p = 2`, as the grid is refined.

## 3. Error norms and order estimation

For each grid level, [`numerics/errors.py`](../src/mms_viscoelastic/numerics/errors.py)
provides:

```text
L2(e)   = sqrt( mean(e_ij^2) )      -- discrete L2 norm
Linf(e) = max(|e_ij|)               -- discrete maximum norm
```

where `e = numerical - exact`.

Given a sequence of mesh sizes `h_k` and corresponding errors `E_k`, the
observed order of accuracy `p` is estimated from the model `E_k ~ C h_k^p`
by a least-squares fit of `log(E_k)` against `log(h_k)`:

```text
log(E_k) = log(C) + p * log(h_k)
```

implemented in `estimate_order` as `numpy.polyfit(log(h), log(E), deg=1)`,
returning the slope `p`.

## 4. `run_convergence_study`

```python
from mms_viscoelastic.numerics import run_convergence_study

report = run_convergence_study([8, 16, 32, 64])
```

For each `n` in `grid_sizes`, the study:

1. builds `Grid2D.from_bounds(nx=n+1, ny=n+1)` on `[0,1]^2`, so `h = dx = dy = 1/n`;
2. evaluates the exact case-12 fields with `evaluate_case_12_exact`;
3. computes `u_num`, `v_num`, `wz_num` from (1)-(3) using the discrete operators;
4. computes `L2` and `Linf` errors of `u_num`, `v_num`, `wz_num` against the
   exact `u`, `v`, `wz`;
5. fits the observed order `p` for each quantity and each norm across all
   requested grids.

The result is a `ConvergenceReport` with one `ConvergenceLevel` per grid and
the fitted orders `l2_orders` / `linf_orders` for `u`, `v`, `wz`.

## 5. Results for `grid_sizes = [8, 16, 32, 64]`

Reproduced with:

```powershell
python examples\09_case12_convergence_verification.py
```

| n  | h        | L2(u)    | Linf(u)  | L2(v)    | Linf(v)  | L2(wz)   | Linf(wz) |
|----|----------|----------|----------|----------|----------|----------|----------|
| 8  | 1.25e-01 | 9.17e-03 | 2.37e-02 | 2.98e-02 | 8.05e-02 | 1.19e-01 | 2.84e-01 |
| 16 | 6.25e-02 | 2.36e-03 | 8.08e-03 | 6.46e-03 | 1.79e-02 | 3.19e-02 | 7.31e-02 |
| 32 | 3.13e-02 | 5.83e-04 | 2.28e-03 | 1.50e-03 | 4.06e-03 | 7.89e-03 | 2.05e-02 |
| 64 | 1.56e-02 | 1.43e-04 | 6.01e-04 | 3.62e-04 | 9.54e-04 | 1.93e-03 | 5.50e-03 |

Fitted orders (log-log fit across the four grids above):

| Quantity | `p` (L2) | `p` (Linf) |
|----------|----------|------------|
| `u`      | 2.0021   | 1.7726     |
| `v`      | 2.1198   | 2.1332     |
| `wz`     | 1.9854   | 1.8902     |

## 6. Interpretation

- All `L2` orders are within `~6%` of the theoretical value `p = 2`,
  confirming second-order accuracy of `derivative_x`/`derivative_y` on the
  case-12 manufactured fields.
- The `Linf` order for `u` (`1.77`) is below `2` because the `Linf` norm is
  dominated by the boundary one-sided stencil, whose error constant is
  larger than the interior central-difference constant
  (`-h^2/3 f'''` vs. `h^2/6 f'''`, see Section 2). This is a well-known
  effect of asymmetric boundary stencils and does not indicate a defect in
  the implementation: the **rate** is still close to `2` and the **L2**
  rate (which averages over the much larger set of interior points) is
  `2.00`.
- The regression test `tests/test_case12_convergence.py` therefore asserts:
  - `L2` order `>= 1.8` for `u`, `v`, `wz` (theoretical `p = 2`, with margin);
  - `Linf` order `>= 1.5` for `u`, `v`, `wz` (looser bound, accounting for
    the boundary-stencil effect above);
  - monotonically decreasing `L2` errors under refinement.

## 7. Reproducing

```powershell
python examples\09_case12_convergence_verification.py
python -m pytest tests\test_case12_convergence.py -v
```
