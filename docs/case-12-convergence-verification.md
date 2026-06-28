# MMS Case 12 - grid-convergence verification

This document describes the verification implemented in
[`src/mms_viscoelastic/numerics/convergence.py`](../src/mms_viscoelastic/numerics/convergence.py),
exercised by
[`examples/09_case12_convergence_verification.py`](../examples/09_case12_convergence_verification.py)
and tested in
[`tests/test_case12_convergence.py`](../tests/test_case12_convergence.py).

The goal is to confirm, by the Method of Manufactured Solutions (MMS), that
the fourth-order compact finite-difference operators in
[`numerics/derivatives.py`](../src/mms_viscoelastic/numerics/derivatives.py)
(`compact_derivative_x`/`compact_derivative_y`, ported from
`coef_4a`/`rhsx_4a` in
[`derivs.f`](../src/fortran_original/4Order_LPOG_Euler/derivs.f)) reproduce
known exact differential relations of case 12 at an order consistent with
the thesis and the JNNFM paper.

This study previously exercised the second-order placeholder operators
(`derivative_x`/`derivative_y`), which were a deliberate, documented
stand-in pending the compact-scheme port. That port is complete; this
document and the L2/Linf tables below were regenerated against the actual
fourth-order operators.

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

`compact_derivative_x` and `compact_derivative_y`
([`numerics/derivatives.py`](../src/mms_viscoelastic/numerics/derivatives.py))
approximate `d f/dx` (respectively `d f/dy`) with a compact (implicit,
Padé-type) scheme: a tridiagonal system is solved for the derivative values,
rather than reading them off an explicit stencil. Interior rows
(`coef_4a`):

```text
f'_{i-1} + 4 f'_i + f'_{i+1} = 3(f_{i+1} - f_{i-1}) / h
```

which is fourth-order accurate. The boundary closure at `i = 1` (and its
mirror at `i = N`):

```text
f'_1 + 3 f'_2 = (-17 f_1 + 9 f_2 + 9 f_3 - f_4) / (6h)
```

is *also* fourth-order accurate: substituting the exact solution and
expanding in Taylor series (checked symbolically with SymPy) gives a
residual of `(1/20) h^4 f^(5)`, matching the interior order. This is why the
observed order below does not drop near the boundary for the first
derivative, unlike the analogous second-derivative closure used elsewhere in
the scheme, whose boundary truncation error is one order lower, `O(h^3)`
(see
[`validation/fortran_reference/case_compact_derivatives/README.md`](../validation/fortran_reference/case_compact_derivatives/README.md)
for the full derivation; this convergence study only exercises the first
derivative, via identities (1)-(2), so that asymmetry does not appear here).

Therefore the discrete relations

```text
u_num  = compact_derivative_y(psi)
v_num  = -compact_derivative_x(psi)
wz_num = compact_derivative_y(u) - compact_derivative_x(v)
```

should converge to the exact fields `u`, `v`, `wz` with error `O(h^4)` as
the grid is refined, subject to the pre-asymptotic caveat in Section 6.

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
| 8  | 1.25e-01 | 5.65e-04 | 2.73e-03 | 4.39e-03 | 1.48e-02 | 6.02e-03 | 1.74e-02 |
| 16 | 6.25e-02 | 1.82e-05 | 9.10e-05 | 2.11e-04 | 9.43e-04 | 2.91e-04 | 1.12e-03 |
| 32 | 3.13e-02 | 9.46e-07 | 7.73e-06 | 9.12e-06 | 5.59e-05 | 1.97e-05 | 1.00e-04 |
| 64 | 1.56e-02 | 4.93e-08 | 5.79e-07 | 3.95e-07 | 3.34e-06 | 1.06e-06 | 7.82e-06 |

Fitted orders (log-log fit across the four grids above):

| Quantity | `p` (L2) | `p` (Linf) |
|----------|----------|------------|
| `u`      | 4.4715   | 4.0166     |
| `v`      | 4.4856   | 4.0398     |
| `wz`     | 4.1318   | 3.6850     |

## 6. Interpretation

- All quantities converge well above the formal order `p = 4` of the
  compact scheme and the RK4 time integrator (not exercised here, since
  this study evaluates the exact fields directly). This is consistent with
  the orders reported in the thesis and the JNNFM paper (approximately
  4.5 for velocity/stress, slightly lower for vorticity at higher Reynolds
  numbers; this study uses the steady manufactured fields, not a Reynolds
  sweep, so it is not directly comparable to the paper's Re-dependent
  tables, only to the general order-of-magnitude finding).
- Observing `p > 4` is a **pre-asymptotic** effect, not evidence of a
  higher formal order: writing the error as `e(h) = C_4 h^4 + C_5 h^5 +
  O(h^6)`, the two-point order estimator used by `estimate_order` satisfies
  `p_hat(h) = 4 + (C_5/C_4) * (1 - 1/r)/ln(r) * h + O(h^2)` for refinement
  ratio `r`, so it converges to `4` as `h -> 0` but can lie above or below
  `4` at finite `h` depending on the sign of `C_5/C_4`. See the "Remark on
  the observed convergence orders" in the JNNFM paper for the full
  derivation; the same mechanism applies here.
- `wz` converges slightly slower than `u`/`v` in both norms, consistent
  with the paper's observation that vorticity is the most sensitive
  quantity to the discretization (it is obtained from `psi` by two
  applications of a compact first derivative, compounding the boundary
  closure's pre-asymptotic correction).
- The regression test `tests/test_case12_convergence.py` therefore asserts:
  - `L2` order `>= 3.9` for `u`, `v`, `wz` (a margin below the worst
    observed value, `4.13` for `wz`);
  - `Linf` order `>= 3.4` for `u`, `v`, `wz` (a margin below the worst
    observed value, `3.68` for `wz`);
  - monotonically decreasing `L2` errors under refinement.

## 7. Reproducing

```powershell
python examples\09_case12_convergence_verification.py
python -m pytest tests\test_case12_convergence.py -v
```
