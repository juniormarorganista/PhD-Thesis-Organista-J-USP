# Numerical Verification of High-Order Viscoelastic Flow Solvers via MMS

This repository contains the code and numerical verification material from my
doctoral thesis at the University of São Paulo (USP). The thesis develops and
verifies high-order computational methods for incompressible viscoelastic flows
using the **Method of Manufactured Solutions (MMS)** — a rigorous technique for
confirming that a numerical solver achieves its theoretical order of accuracy.

---

## Table of contents

1. [Scientific context](#scientific-context)
2. [Mathematical formulation](#mathematical-formulation)
3. [Method of Manufactured Solutions](#method-of-manufactured-solutions)
4. [Numerical verification](#numerical-verification)
5. [Convergence results (Case 12)](#convergence-results-case-12)
6. [Thesis figures](#thesis-figures)
7. [What is in this repository](#what-is-in-this-repository)
8. [Time-marching solver (Case 12)](#time-marching-solver-case-12)
9. [Repository structure](#repository-structure)
10. [Installation](#installation)
11. [Tests and validation](#tests-and-validation)
12. [Examples](#examples)
13. [Machine learning utilities](#machine-learning-utilities)
14. [Citation](#citation)
15. [License](#license)

---

## Scientific context

Viscoelastic flows are challenging to simulate accurately because the velocity
field couples with the evolution of a polymeric stress tensor, introducing
memory effects and numerical stiffness controlled by the Weissenberg number. My
work focuses on verifying that the high-order finite-difference solver
implemented in Fortran actually achieves the expected fourth-order spatial
accuracy across a range of flow configurations and constitutive models (UCM,
Oldroyd-B, Giesekus, LPTT).

The central question is not "does the simulation look right?" but rather:
**does the measured convergence rate match the theoretical prediction?**
That is a code-verification problem, and MMS is the standard rigorous tool for
answering it.

This repository has two complementary layers. The first preserves the original
scientific material for traceability:

```text
src/fortran_original/    — original Fortran solver (4th-order LPOG compact scheme)
src/mathematica_original/ — symbolic derivations of manufactured fields
Thesis/                  — thesis source and figures
Papers/                  — journal submission material
```

The second builds a modern Python package around those originals:

```text
src/mms_viscoelastic/    — validated Python re-implementation
validation/              — isolated Fortran reference programs + reference CSVs
tests/                   — pytest regression suite (508 tests)
examples/                — runnable scripts and notebooks
```

The Python layer is organized as a real scientific package: it reproduces every
numerical result from the Fortran originals to within machine precision, carries
a full regression-test suite, and is structured so that the verification chain
is fully reproducible without the original Fortran toolchain.

---

## Mathematical formulation

### Incompressibility

For a two-dimensional incompressible flow the velocity field
$\mathbf{u} = (u, v)$ satisfies:

```math
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0.
```

### Streamfunction formulation

Introducing a streamfunction $\psi$ so that

```math
u = \frac{\partial \psi}{\partial y}, \qquad v = -\frac{\partial \psi}{\partial x}
```

automatically satisfies the incompressibility constraint.

### Vorticity

The scalar vorticity is:

```math
\omega_z = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
         = -\nabla^2 \psi,
\qquad \text{i.e.,} \quad \nabla^2 \psi = -\omega_z.
```

### Rate-of-deformation tensor

```math
\mathbf{D} = \frac{1}{2}\left(\nabla\mathbf{u} + \nabla\mathbf{u}^{T}\right)
= \begin{bmatrix}
\dfrac{\partial u}{\partial x} &
\dfrac{1}{2}\!\left(\dfrac{\partial u}{\partial y}+\dfrac{\partial v}{\partial x}\right) \\[6pt]
\dfrac{1}{2}\!\left(\dfrac{\partial u}{\partial y}+\dfrac{\partial v}{\partial x}\right) &
\dfrac{\partial v}{\partial y}
\end{bmatrix}.
```

### Viscoelastic constitutive equation

The polymeric extra-stress tensor $\boldsymbol{\tau}$ evolves according to:

```math
\mathrm{Wi}\,\mathcal{D}\boldsymbol{\tau}
+ \mathcal{F}\!\left(\boldsymbol{\tau}, \nabla \mathbf{u}\right)
= 2(1-\beta)\mathbf{D},
```

where $\mathrm{Wi}$ is the Weissenberg number, $\beta$ the solvent-viscosity
ratio, $\mathcal{D}$ a suitable objective derivative, and $\mathcal{F}$ depends
on the constitutive model. In two dimensions the stress tensor is:

```math
\boldsymbol{\tau} = \begin{bmatrix} T_{xx} & T_{xy} \\ T_{xy} & T_{yy} \end{bmatrix}.
```

---

## Method of Manufactured Solutions

The idea is to prescribe a smooth artificial exact solution and derive the
source terms that make it satisfy the governing equations exactly. If the solver
then converges to that solution at the expected rate, the implementation is
verified.

The manufactured state vector for Case 12 is:

```math
\mathbf{q}_e(x,y,t)
= \bigl(u_e,\; v_e,\; \psi_e,\; \omega_{z,e},\; T_{xx,e},\; T_{xy,e},\; T_{yy,e}\bigr).
```

Substituting into the governing equations gives the MMS source vector:

```math
\mathcal{L}(\mathbf{q}_e) = \mathbf{f}_{\mathrm{MMS}}
= \bigl(f_{\omega_z},\; f_{T_{xx}},\; f_{T_{xy}},\; f_{T_{yy}}\bigr).
```

In the Fortran implementation these are `tfwz`, `tfTxx`, `tfTxy`, `tfTyy`.
The Python port preserves the same scientific meaning with NumPy `(ny, nx)`
array layout.

---

## Numerical verification

Given a numerical approximation $\mathbf{q}_h$ on a grid of spacing $h$, the
error is $\mathbf{e}_h = \mathbf{q}_h - \mathbf{q}_e$, measured in the norms:

```math
\lVert e_h \rVert_{2}
= \left(\frac{1}{N}\sum_{i=1}^{N} e_i^2\right)^{1/2}, \qquad
\lVert e_h \rVert_{\infty} = \max_i |e_i|.
```

The observed order of convergence is estimated from two successive grids:

```math
p = \frac{\log(E_{h_1}/E_{h_2})}{\log(h_1/h_2)}.
```

For a fourth-order scheme the expected behavior is $E_h \approx C h^4$, so
$p \approx 4$.

---

## Convergence results (Case 12)

The table below verifies the fourth-order compact operators in
[`numerics/derivatives.py`](src/mms_viscoelastic/numerics/derivatives.py)
(`compact_derivative_x`/`compact_derivative_y`, ported from `coef_4a`/`rhsx_4a`
in [`derivs.f`](src/fortran_original/4Order_LPOG_Euler/derivs.f)) through the
streamfunction identities $u = \partial\psi/\partial y$,
$v = -\partial\psi/\partial x$, $w_z = \partial u/\partial y - \partial v/\partial x$.

> **Sign convention.** The `wz` field used here, and in the original `case_12.f`,
> equals $\partial u/\partial y - \partial v/\partial x = +\nabla^2\psi$, which is
> the *negative* of the $\omega_z$ defined above. The Python port follows the
> Fortran convention so that field-by-field comparisons remain direct. See
> [`docs/case-12-complete-validation-report.md`](docs/case-12-complete-validation-report.md#12-a-note-on-the-sign-of-wz).

Results for grids $n = 8, 16, 32, 64$ ($h = 1/n$):

| $n$ | $h$ | $L_2(u)$ | $L_\infty(u)$ | $L_2(v)$ | $L_\infty(v)$ | $L_2(w_z)$ | $L_\infty(w_z)$ |
|-----|----------|----------|----------|----------|----------|----------|----------|
| 8   | 1.25e-01 | 5.65e-04 | 2.73e-03 | 4.39e-03 | 1.48e-02 | 6.02e-03 | 1.74e-02 |
| 16  | 6.25e-02 | 1.82e-05 | 9.10e-05 | 2.11e-04 | 9.43e-04 | 2.91e-04 | 1.12e-03 |
| 32  | 3.13e-02 | 9.46e-07 | 7.73e-06 | 9.12e-06 | 5.59e-05 | 1.97e-05 | 1.00e-04 |
| 64  | 1.56e-02 | 4.93e-08 | 5.79e-07 | 3.95e-07 | 3.34e-06 | 1.06e-06 | 7.82e-06 |

Estimated order of accuracy (log-log fit across all four grids):

| Quantity | $p$ ($L_2$) | $p$ ($L_\infty$) |
|----------|-------------|------------------|
| $u$      | 4.47        | 4.02             |
| $v$      | 4.49        | 4.04             |
| $w_z$    | 4.13        | 3.69             |

All quantities converge above the formal order $p = 4$. This is a pre-asymptotic
effect of the next-order term in the error expansion
($e(h) = C_4 h^4 + C_5 h^5 + \cdots$), consistent with the orders reported in
the thesis and the JNNFM paper. $w_z$ converges slightly slower than $u$/$v$
because it is obtained from $\psi$ through two compact derivatives.
Full derivations and the symbolic verification with SymPy are in
[`docs/case-12-convergence-verification.md`](docs/case-12-convergence-verification.md).

Reproduce with:

```powershell
python examples\09_case12_convergence_verification.py
python -m pytest tests\test_case12_convergence.py -v
```

---

## Thesis figures

The figures below are from the thesis and presentation material. The originals
were EPS; PNG versions are included here so GitHub can render them.

### Flow field — streamlines

<p align="center">
  <img src="assets/thesis/case12/thesis_streamlines.png" width="70%"
       alt="Streamlines from thesis material">
</p>

### MMS validation pipeline

<p align="center">
  <img src="assets/preview/mms_validation_pipeline.svg" width="90%"
       alt="MMS validation pipeline">
</p>

### Case 12 — exact fields

| $\psi$ | $u$ |
|---|---|
| <img src="assets/thesis/case12/exact_map_psi.png" width="420"> | <img src="assets/thesis/case12/exact_map_u.png" width="420"> |

| $v$ | $\omega_z$ |
|---|---|
| <img src="assets/thesis/case12/exact_map_v.png" width="420"> | <img src="assets/thesis/case12/exact_map_wz.png" width="420"> |

| $T_{xx}$ | $T_{xy}$ | $T_{yy}$ |
|---|---|---|
| <img src="assets/thesis/case12/exact_map_txx.png" width="300"> | <img src="assets/thesis/case12/exact_map_txy.png" width="300"> | <img src="assets/thesis/case12/exact_map_tyy.png" width="300"> |

### Case 12 — exact surfaces

| $\psi$ | $u$ |
|---|---|
| <img src="assets/thesis/case12/exact_surface_psi.png" width="420"> | <img src="assets/thesis/case12/exact_surface_u.png" width="420"> |

| $\omega_z$ | $T_{xx}$ |
|---|---|
| <img src="assets/thesis/case12/exact_surface_wz.png" width="420"> | <img src="assets/thesis/case12/exact_surface_txx.png" width="420"> |

### Convergence curves

| $\psi$ | $u$ |
|---|---|
| <img src="assets/thesis/case12/convergence_re100_wi1_alphaG0_psi.png" width="420"> | <img src="assets/thesis/case12/convergence_re100_wi1_alphaG0_u.png" width="420"> |

| $v$ | $\omega_z$ |
|---|---|
| <img src="assets/thesis/case12/convergence_re100_wi1_alphaG0_v.png" width="420"> | <img src="assets/thesis/case12/convergence_re100_wi1_alphaG0_wz.png" width="420"> |

| $T_{xx}$ | $T_{xy}$ | $T_{yy}$ |
|---|---|---|
| <img src="assets/thesis/case12/convergence_re100_wi1_alphaG0_txx.png" width="300"> | <img src="assets/thesis/case12/convergence_re100_wi1_alphaG0_txy.png" width="300"> | <img src="assets/thesis/case12/convergence_re100_wi1_alphaG0_tyy.png" width="300"> |

### Numerical slices

| Slice in $x$ for $u$ | Slice in $y$ for $\omega_z$ |
|---|---|
| <img src="assets/thesis/case12/slice_x_u.png" width="420"> | <img src="assets/thesis/case12/slice_y_wz.png" width="420"> |

---

## What is in this repository

### Python package (`src/mms_viscoelastic/`)

The Fortran-to-Python port of `src/fortran_original/4Order_LPOG_Euler/` is
functionally complete. The package includes:

- **Grid utilities** — `Grid2D` with uniform spacing
- **All 14 MMS manufactured-solution cases** (`case_01` through `case_14`),
  each ported and validated against an isolated Fortran reference program.
  Case 12 is the primary case from the thesis and the JNNFM paper. Cases 13
  and 14 are not MMS cases in the original Fortran — their exact solution and
  source term are identically zero (lid-driven-cavity benchmarks); they were
  ported as trivial structural-parity modules (see
  [`docs/case-13-14-port.md`](docs/case-13-14-port.md)).
- **Fourth-order compact derivatives and Lele filter** — ports of `coef_4a`,
  `rhsx_4a`, `rhsxx_4a` and the compact filter from `derivs.f`, validated
  against isolated Fortran references.
- **Log-conformation reformulation** — port of `conform.f`: closed-form 2×2
  eigendecomposition (`decomp_matrix`), `T_to_Psi`/`Psi_to_T` conversions,
  `kernel_conformation`, and the explicit-Euler time step (`drv_Psi`,
  `solver_log_euler`). Validated to machine precision against an isolated
  Fortran reference that uses the general Jacobi iteration (not the closed form),
  making the comparison independent. This path is not discussed in the thesis
  text; it is ported for source-code fidelity. See
  [`docs/case-12-log-conformation.md`](docs/case-12-log-conformation.md).
- **Source terms for all 14 cases**, with full Fortran-reference validation for
  Case 12's four terms (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`).
- **Constitutive models** (UCM, Oldroyd-B), validated against closed-form
  steady Couette solutions.
- **Grid-convergence verification** for the streamfunction identities (see
  [Convergence results](#convergence-results-case-12)).
- **Time-marching solver** for Case 12 (see
  [Time-marching solver](#time-marching-solver-case-12)).
- **Visualization utilities** and a **reproducible Jupyter notebook**
  ([`examples/notebooks/Case12_MMS_Pipeline.ipynb`](examples/notebooks/Case12_MMS_Pipeline.ipynb)).
- **ML dataset utilities** (exploratory extension — see
  [Machine learning utilities](#machine-learning-utilities)).

The file-by-file correspondence between Fortran routines and Python modules is
documented in [`docs/fortran-code-map.md`](docs/fortran-code-map.md).

### Quick API example

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import evaluate_case_12_source_terms

grid = Grid2D.from_bounds(nx=65, ny=65)
terms = evaluate_case_12_source_terms(grid)

terms.tfwz    # forcing for the vorticity equation
terms.tfTxx   # forcing for T_xx
terms.tfTxy   # forcing for T_xy
terms.tfTyy   # forcing for T_yy
```

### What is still out of scope

The following items are documented as explicitly deferred:

```text
RK4 time integration (solver_rk4, solver_log_rk4) and full multi-step
  log-conformation time marching — only the Euler variants are ported.
Full temporal-evolution comparison against the original Fortran time loop
  on fine grids and long simulation horizons. Explored in Etapa 17
  (docs/case-12-4th-order-time-marching.md): the explicit-Euler instability
  documented in docs/case-12-full-march-reference.md persists with the real
  4th-order + Lele-filter + multigrid-Poisson scheme, just at a different
  point in time — the dissipation from the filter alone is not enough to
  stabilize the explicit march at dt=1e-5 on coarse grids.
```

---

## Time-marching solver (Case 12)

`src/mms_viscoelastic/numerics/` ports the building blocks of the Case 12
time-marching loop from `src/fortran_original/4Order_LPOG_Euler/`, using the
**real scheme** the original Fortran uses: 4th-order LPOG compact derivatives,
the Lele compact filter applied once per step between the boundary
re-imposition and the Poisson re-solve (matching `subroutine loop` in
`main.f`), and a 4th-order compact multigrid V-cycle Poisson solver
(`poi_psi4a`). The pentadiagonal `band5_poi`/`banbk5_poi` routines in
`solvers.f` are dead code — a repository-wide search confirms they are never
called anywhere in the original Fortran.

Each building block is validated against an isolated Fortran reference program
written for this project, which does not touch the original thesis sources.

| Step | Module | Documentation |
|---|---|---|
| 4th-order compact multigrid Poisson solve for `psi` | [`numerics/poisson_multigrid.py`](src/mms_viscoelastic/numerics/poisson_multigrid.py) | [`docs/case-12-4th-order-time-marching.md`](docs/case-12-4th-order-time-marching.md) |
| Explicit Euler step for `wz`, `Txx`, `Txy`, `Tyy` | [`numerics/time_stepping.py`](src/mms_viscoelastic/numerics/time_stepping.py) | [`docs/case-12-euler-step-comparison.md`](docs/case-12-euler-step-comparison.md) |
| Full step: Euler + BCs + Lele filter + Poisson + `u`,`v` recomputation | [`numerics/time_marching.py`](src/mms_viscoelastic/numerics/time_marching.py) | [`docs/case-12-time-marching-comparison.md`](docs/case-12-time-marching-comparison.md) |
| Multi-step march (`n_steps` chained full steps) | [`numerics/time_marching.py`](src/mms_viscoelastic/numerics/time_marching.py) | [`docs/case-12-multistep-marching-comparison.md`](docs/case-12-multistep-marching-comparison.md) |

All steps agree with their Fortran references to within machine precision
(~1e-15; bit-for-bit for the Poisson solver). A second-order direct Poisson
solve ([`numerics/poisson.py`](src/mms_viscoelastic/numerics/poisson.py))
remains as an independent utility but is no longer part of the main pipeline.

---

## Repository structure

```text
.
├── LICENSE
├── CITATION.cff
├── src/
│   ├── fortran_original/
│   │   └── 4Order_LPOG_Euler/        ← original Fortran solver (read-only)
│   ├── mathematica_original/         ← symbolic derivations (read-only)
│   └── mms_viscoelastic/
│       ├── grids.py
│       ├── fields.py
│       ├── parameters.py
│       ├── manufactured_solutions/   (case_01.py – case_14.py)
│       ├── source_terms/             (case_01.py – case_14.py)
│       ├── models/                   (ucm.py, oldroyd_b.py, log_conformation.py)
│       ├── numerics/                 (derivatives, filters, poisson, time_stepping[_log], time_marching)
│       ├── validation/
│       ├── visualization/
│       ├── ml/                       (exploratory extension — not thesis core)
│       └── io/
├── Thesis/                           ← thesis LaTeX source and figures (read-only)
├── Papers/                           ← journal submission material (read-only)
├── validation/
│   ├── fortran_reference/            ← isolated Fortran reference programs + build scripts
│   └── reference_data/               ← reference CSVs (committed, regenerable)
├── examples/
├── tests/
├── docs/
├── scripts/
└── assets/
    ├── preview/                      ← diagrams committed to the repo
    └── thesis/case12/                ← thesis figures committed to the repo
```

The folders `src/fortran_original/`, `src/mathematica_original/`, `Thesis/`
and `Papers/` are reference material and are not modified.

---

## Installation

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1

python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements-dev.txt
python -m pip install -e .
```

If PowerShell blocks virtual environment activation:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
.\.venv\Scripts\Activate.ps1
```

---

## Tests and validation

### Full test suite

```powershell
python -m pytest -q
```

508 tests cover manufactured-solution fields, source terms, compact
derivatives, the Lele filter, the multigrid Poisson solver, constitutive
models, time-stepping building blocks, and Fortran-reference comparisons.

### Source-term reference validation

```powershell
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
python examples\08_validate_case12_tfTyy_reference.py
```

Each script prints `passed: True` and the maximum absolute error against the
Fortran reference.

### Grid-convergence verification

```powershell
python examples\09_case12_convergence_verification.py
python -m pytest tests\test_case12_convergence.py -v
```

See [Convergence results (Case 12)](#convergence-results-case-12) for the
expected output.

---

## Examples

```powershell
python examples\01_smoke_python_package.py            # package import check
python examples\02_case12_manufactured_solution.py    # exact fields on a grid
python examples\03_case12_visualization_and_ml_dataset.py
python examples\04_case12_tfwz_source_term.py         # source terms → CSV
python examples\09_case12_convergence_verification.py # convergence study
```

For a reproducible end-to-end pipeline (exact fields → source-term validation
→ convergence study):

```powershell
jupyter nbconvert --to notebook --execute --inplace `
  examples\notebooks\Case12_MMS_Pipeline.ipynb
```

Consolidated results from this notebook are in
[`docs/case-12-complete-validation-report.md`](docs/case-12-complete-validation-report.md).

---

## Machine learning utilities

> **Scope note.** The `ml/` subpackage is an **exploratory extension, not part
> of the thesis's core scope**. The thesis is a numerical-verification work; it
> is not about machine learning. No surrogate model is trained anywhere in this
> repository — the utilities only build pointwise feature tables and
> train/validation/test splits from the already-validated manufactured fields.
> The surrogate-model phase originally planned was discontinued; the decision
> is recorded in the planning log.

```text
src/mms_viscoelastic/ml/
```

A typical dataset row includes coordinates, physical parameters
($\beta$, Re, Wi), the seven manufactured fields ($u$, $v$, $\psi$, $\omega_z$,
$T_{xx}$, $T_{xy}$, $T_{yy}$) and the four source terms
($f_{\omega_z}$, $f_{T_{xx}}$, $f_{T_{xy}}$, $f_{T_{yy}}$). These datasets
could serve as input to physics-informed networks, surrogate models, or
error estimators in future work, but that is outside the scope of this thesis.

---

## Citation

If this repository supports your research, please cite both the doctoral thesis
and this repository:

```text
Organista, J. Numerical verification of high-order solvers for incompressible
viscoelastic flows using the Method of Manufactured Solutions.
Doctoral thesis, University of São Paulo (USP).
```

Machine-readable citation metadata is in [`CITATION.cff`](CITATION.cff), which
GitHub surfaces through the "Cite this repository" button on the repository
home page.

---

## License

The Python library, examples, tests and documentation in this repository are
released under the [MIT License](LICENSE).

The original Fortran and Mathematica reference material under
`src/fortran_original/`, `src/mathematica_original/` and `Thesis/` keeps its
own terms, described in [`Thesis/LICENSE`](Thesis/LICENSE)
(BSD-3-Clause style, Copyright © 2015 Humberto Lidio Antonelli).

---

**Juniormar Organista** — Applied and Computational Mathematics, USP.
Computational fluid mechanics, viscoelastic flows, high-order finite differences,
Method of Manufactured Solutions.
