# Case 12 - Complete MMS Validation Report

This report consolidates everything that has been verified so far for **Case
12** of the original Fortran solver
([`src/fortran_original/4Order_LPOG_Euler/case_12.f`](../src/fortran_original/4Order_LPOG_Euler/case_12.f))
and its Python port in `mms_viscoelastic`. Case 12 corresponds to the
manufactured-solution test labeled `MMS_12` in the thesis figures (see
[`Thesis/figures/Case12/`](../Thesis/figures/Case12/)), and it is the case
chosen as the entry point for the Python modernization because it exercises
the full state vector - velocity, streamfunction, vorticity and the three
extra-stress components - with a single closed-form manufactured solution.

The report brings together three independent pieces of evidence:

1. the exact manufactured fields and the algebraic identities they satisfy;
2. the Fortran-reference validation of the four manufactured source terms
   (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`);
3. the grid-convergence study of the finite-difference operators applied to
   the stream-function identities.

Together, they show that the Python port reproduces the Fortran
manufactured solution exactly (to floating-point precision) and that the
discrete operators used to check it converge at the expected order. The
reproducible notebook
[`examples/notebooks/Case12_MMS_Pipeline.ipynb`](../examples/notebooks/Case12_MMS_Pipeline.ipynb)
runs all three checks end to end and produces the figures referenced below.

---

## 1. The exact manufactured solution

The exact fields are implemented in
[`evaluate_case_12_exact`](../src/mms_viscoelastic/manufactured_solutions/case_12.py),
a direct port of the Fortran subroutine `exact_solutions_case_12(t)`. On the
unit square $[0,1]^2$, with $x_p(x) = (x-1)x$ and decay factor
$\mathrm{dec} = e^{a_t t}$, the fields are:

```math
\psi
=
\frac{x_p(x)\, y \, \sin(\pi x) \sin(\pi y)}{\mathrm{dec}},
```

```math
u
=
\frac{x_p(x) \sin(\pi x)\big(\pi y \cos(\pi y) + \sin(\pi y)\big)}{\mathrm{dec}},
```

```math
v
=
-\,
\frac{y \big(\pi x_p(x) \cos(\pi x) + (2x-1)\sin(\pi x)\big)\sin(\pi y)}{\mathrm{dec}},
```

```math
w_z
=
\frac{2\Big(\pi x_p(x)\cos(\pi y)\sin(\pi x)
+ y\big(\pi(2x-1)\cos(\pi x) + (1-\pi^2 x_p(x))\sin(\pi x)\big)\sin(\pi y)\Big)}{\mathrm{dec}}.
```

The three stress components share the same separable structure,

```math
T_{xx}, T_{xy}, T_{yy}
=
\frac{(1-\beta)\, S_x(x)\, S_y(y)}{\mathrm{dec}},
\qquad
S_x(x) = a_1 + a_3\cos(\pi x) + a_2\sin(\pi x),
\quad
S_y(y) = b_1 + b_3\cos(\pi y) + b_2\sin(\pi y),
```

with independent coefficient sets $(a_1,a_2,a_3,b_1,b_2,b_3)$ for each of
`txx`, `txy`, `tyy` (all equal to one by default, matching `par.nn` in the
original code).

At $t=0$ on a $65\times 65$ grid these fields take the ranges shown in the
notebook:

| field | min | max |
|---|---|---|
| `u`   | $-3.4775\times 10^{-1}$ | $7.8540\times 10^{-1}$ |
| `v`   | $-4.4851\times 10^{-1}$ | $4.4851\times 10^{-1}$ |
| `psi` | $-1.4478\times 10^{-1}$ | $0.0$ |
| `wz`  | $-3.6387$ | $4.8522$ |
| `txx`, `txy`, `tyy` | $\approx 0$ | $2.9142$ |

### 1.1 Stream-function identities

These five fields are not independent. Differentiating the closed-form
expressions above and simplifying with SymPy shows, term by term, that the
following three relations hold **identically** (the symbolic residual is
exactly `0`, for arbitrary $x$, $y$, $t$ and $a_t$):

```math
u - \frac{\partial \psi}{\partial y} = 0,
\qquad
v + \frac{\partial \psi}{\partial x} = 0,
\qquad
w_z - \left(\frac{\partial u}{\partial y} - \frac{\partial v}{\partial x}\right) = 0.
```

The third relation is equivalent to $w_z = \nabla^2\psi$, since
$\partial u/\partial y - \partial v/\partial x = \partial^2\psi/\partial y^2 +
\partial^2\psi/\partial x^2$ once the first two relations are substituted.

These identities have no modeling content - they are algebraic consequences
of the chosen $\psi$, $u$, $v$ and $w_z$ - but they are exactly the relations
the discrete operators must reproduce, which is what Section 3 checks.

### 1.2 A note on the sign of `wz`

The README's [Mathematical formulation](../README.md#mathematical-formulation)
section defines vorticity through the standard planar curl,

```math
\omega_z = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
= -\nabla^2\psi.
```

Case 12's `wz`, as written in `case_12.f` and reproduced above, satisfies the
**opposite** relation, $w_z = \partial u/\partial y - \partial v/\partial x =
+\nabla^2\psi = -\omega_z$. Both sign conventions for "vorticity" are common
in the viscoelastic-flow literature, and the Fortran code consistently uses
$w_z = +\nabla^2\psi$ throughout (the manufactured field, the time-stepping
unknown, and the source term `tfwz` are all written in this convention).

The Python port preserves this convention for one reason: fidelity to the
original solver. Re-deriving Case 12 with the opposite sign would require
renegotiating every formula in `case_12.f`, including the source terms, and
would break the direct field-by-field comparison against the Fortran
reference that is the whole point of this verification effort. The price of
that choice is that `wz` in this package is the negative of the $\omega_z$ in
the README's general formulation - a fact now documented here, in the
[convergence-verification note](case-12-convergence-verification.md#1-exact-differential-identities-satisfied-by-case-12),
and in the README's
[Convergence results](../README.md#convergence-results-case-12) section, so
that nobody mixes the two conventions when extending the package to other
cases.

---

## 2. Source-term validation against the Fortran reference

The manufactured source terms `tfwz`, `tfTxx`, `tfTxy`, `tfTyy` are the
right-hand sides that make the exact fields of Section 1 satisfy the
governing PDEs exactly. They are implemented in
[`evaluate_case_12_source_terms`](../src/mms_viscoelastic/source_terms/case_12.py)
and were ported term by term from `source_term_case_12(t)`.

Each term is checked against an **independent** implementation: a small
Fortran program
([`validation/fortran_reference/case12_<term>/`](../validation/fortran_reference/))
evaluates the same closed-form expression on a $3\times3$ (and, for some
terms, $5\times5$) grid and writes the result to a CSV file
([`validation/reference_data/`](../validation/reference_data/)). The Python
values are then compared point by point with
[`compare_case12_*_reference`](../src/mms_viscoelastic/validation/source_term_reference.py).

Running all four checks (`examples/05` through `examples/08`, also reproduced
in the notebook) on the $3\times3$ default grid gives:

| term | `max_abs_error` | `max_rel_error` | acceptance | result |
|---|---|---|---|---|
| `tfwz`  | $7.105427\times 10^{-15}$ | $1.439861\times 10^{-16}$ | $< 10^{-11}$ | **passed** |
| `tfTxx` | $4.440892\times 10^{-16}$ | $2.066104\times 10^{-16}$ | $< 10^{-10}$ | **passed** |
| `tfTxy` | $8.881784\times 10^{-16}$ | $1.498532\times 10^{-16}$ | $< 10^{-10}$ | **passed** |
| `tfTyy` | $1.732348\times 10^{-20}$ | $1.732348\times 10^{-20}$ | $< 10^{-10}$ | **passed** |

All four errors sit at the level of double-precision round-off
($10^{-15}$-$10^{-20}$), several orders of magnitude below the acceptance
thresholds documented in
[`docs/case-12-tfwz-fortran-reference-validation.md`](case-12-tfwz-fortran-reference-validation.md),
[`docs/case-12-tfTxx-fortran-reference-validation.md`](case-12-tfTxx-fortran-reference-validation.md),
[`docs/case-12-tfTxy-fortran-reference-validation.md`](case-12-tfTxy-fortran-reference-validation.md)
and
[`docs/case-12-tfTyy-fortran-reference-validation.md`](case-12-tfTyy-fortran-reference-validation.md).
The differences are consistent with accumulated floating-point rounding
between two independently typed implementations of the same algebraic
expression, not with a modeling discrepancy. With `tfTyy` now validated, all
four Case 12 source terms have been verified against the Fortran reference -
the Python port reproduces `source_term_case_12(t)` exactly.

---

## 3. Grid-convergence of the discrete stream-function identities

The third check, implemented in
[`run_convergence_study`](../src/mms_viscoelastic/numerics/convergence.py)
and documented in full in
[`docs/case-12-convergence-verification.md`](case-12-convergence-verification.md),
takes the identities of Section 1.1 and replaces the analytical derivatives
with the second-order operators `derivative_x` / `derivative_y` from
[`numerics/derivatives.py`](../src/mms_viscoelastic/numerics/derivatives.py):

```text
u_num  = derivative_y(psi)
v_num  = -derivative_x(psi)
wz_num = derivative_y(u) - derivative_x(v)
```

Because Section 1.1 established that the corresponding *exact* relations have
zero residual, any discrepancy between `u_num`/`v_num`/`wz_num` and the exact
`u`/`v`/`wz` is pure discretization error. For grids with $n = 8, 16, 32, 64$
cells per direction ($h = 1/n$):

| $n$ | $h$ | $L_2(u)$ | $L_\infty(u)$ | $L_2(v)$ | $L_\infty(v)$ | $L_2(w_z)$ | $L_\infty(w_z)$ |
|----|----------|----------|----------|----------|----------|----------|----------|
| 8  | $1.25\times10^{-1}$ | $9.17\times10^{-3}$ | $2.37\times10^{-2}$ | $2.98\times10^{-2}$ | $8.05\times10^{-2}$ | $1.19\times10^{-1}$ | $2.84\times10^{-1}$ |
| 16 | $6.25\times10^{-2}$ | $2.36\times10^{-3}$ | $8.08\times10^{-3}$ | $6.46\times10^{-3}$ | $1.79\times10^{-2}$ | $3.19\times10^{-2}$ | $7.31\times10^{-2}$ |
| 32 | $3.13\times10^{-2}$ | $5.83\times10^{-4}$ | $2.28\times10^{-3}$ | $1.50\times10^{-3}$ | $4.06\times10^{-3}$ | $7.89\times10^{-3}$ | $2.05\times10^{-2}$ |
| 64 | $1.56\times10^{-2}$ | $1.43\times10^{-4}$ | $6.01\times10^{-4}$ | $3.62\times10^{-4}$ | $9.54\times10^{-4}$ | $1.93\times10^{-3}$ | $5.50\times10^{-3}$ |

The observed order of accuracy $p$, estimated by a log-log least-squares fit
across all four grids ($\log E_h = \log C + p \log h$), is:

| quantity | $p_{L_2}$ | $p_{L_\infty}$ |
|---|---|---|
| $u$  | 2.0021 | 1.7726 |
| $v$  | 2.1198 | 2.1332 |
| $w_z$ | 1.9854 | 1.8902 |

All three $L_2$ orders sit within 6% of the theoretical value $p=2$ expected
for the central-difference stencil used in the interior. The $L_\infty$ order
for $u$ (1.77) is the only value noticeably below 2; this is explained by the
boundary one-sided stencil, whose truncation-error constant is larger than
the interior central-difference constant ($h^2/3$ vs. $h^2/6$, see
[Section 2 of the convergence note](case-12-convergence-verification.md#2-discrete-operators-and-their-truncation-error)).
It is an expected feature of asymmetric boundary stencils, not a defect: the
$L_2$ norm, which is dominated by the much larger set of interior points,
still recovers $p\approx2.00$.

The regression tests in
[`tests/test_case12_convergence.py`](../tests/test_case12_convergence.py)
encode these observations as `L2 order >= 1.8` and `Linf order >= 1.5` for
`u`, `v`, `wz`, plus a monotonic-decrease check on the $L_2$ errors.

---

## 4. Reproducing this report

All three checks, plus the visualizations referenced above, are reproduced in
a single notebook:

```powershell
jupyter nbconvert --to notebook --execute --inplace `
  examples\notebooks\Case12_MMS_Pipeline.ipynb
```

Equivalently, each check can be run on its own:

```powershell
# Section 1 - exact fields
python examples\02_case12_manufactured_solution.py

# Section 2 - source-term validation against Fortran
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
python examples\08_validate_case12_tfTyy_reference.py

# Section 3 - grid-convergence study
python examples\09_case12_convergence_verification.py

# Full test suite (41 tests at the time of writing)
python -m pytest -q
```

---

## 5. Summary

- The exact fields of Case 12 satisfy the stream-function identities
  $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$,
  $w_z=\partial u/\partial y-\partial v/\partial x=\nabla^2\psi$ exactly
  (verified symbolically).
- All four manufactured source terms (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`) match
  the independent Fortran reference to within double-precision round-off,
  well inside the documented acceptance thresholds.
- The second-order finite-difference operators reproduce the stream-function
  identities with $L_2$ order $\approx 2.0$ across grids $8\to16\to32\to64$.
- The sign of `wz` in this package is $-\omega_z$ relative to the general
  vorticity definition in the README; this is a deliberate choice to keep the
  Python port aligned with `case_12.f`, and is documented here and in the
  convergence note to avoid confusion in future cases.

Taken together, these three results are the current body of MMS verification
evidence for Case 12 in the Python port.
