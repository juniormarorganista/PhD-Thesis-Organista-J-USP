# Thesis USP Codes - MMS Verification and Python Modernization

Repository associated with the doctoral thesis of **Juniormar Organista**, focused on the numerical verification of high-order solvers for incompressible viscoelastic flows using the **Method of Manufactured Solutions (MMS)**.

---

## Table of contents

1. [Project vision](#project-vision)
2. [Scientific context](#scientific-context)
3. [Mathematical formulation](#mathematical-formulation)
4. [Method of Manufactured Solutions](#method-of-manufactured-solutions)
5. [Numerical verification](#numerical-verification)
6. [Thesis figures in this README](#thesis-figures-in-this-readme)
7. [Current Python status](#current-python-status)
8. [Case 12 status](#case-12-status)
9. [Repository structure](#repository-structure)
10. [Installation](#installation)
11. [Tests and validation](#tests-and-validation)
12. [Examples](#examples)
13. [Generated outputs](#generated-outputs)
14. [Machine Learning preparation](#machine-learning-preparation)
15. [Local Git workflow](#local-git-workflow)
16. [GitHub About](#github-about)
17. [Citation](#citation)
18. [Author](#author)
19. [README maintenance rule](#readme-maintenance-rule)

---

## Project vision

The project has two complementary goals.

First, it preserves the scientific material associated with the doctoral thesis. The original sources remain available for traceability, historical inspection and comparison:

```text
src/fortran_original/
src/mathematica_original/
thesis/original/
papers/original/
```

Second, it builds a modern Python package:

```text
src/mms_viscoelastic/
```

The Python layer is being organized as a real scientific package with tests, validation data, examples, source-term ports, visualization utilities and future ML-ready datasets.

The long-term objective is to transform a thesis codebase into a reproducible scientific computing framework.

---

## Scientific context

The thesis concerns the numerical verification of high-order computational solvers for **incompressible viscoelastic flows**.

Viscoelastic flows are challenging because the flow field is coupled with polymeric stress evolution. In practice, the solver must handle:

- velocity components;
- streamfunction;
- scalar vorticity;
- stress tensor components;
- nonlinear constitutive terms;
- stiffness associated with the Weissenberg number;
- high-order spatial discretization;
- manufactured forcing terms.

The central scientific question is:

> Given an exact manufactured solution, does the implemented numerical method recover the expected theoretical order of accuracy?

This is a code-verification problem, not merely a simulation problem.

---

## Mathematical formulation

### Incompressibility

For an incompressible flow, the velocity field satisfies:

```math
\nabla \cdot \mathbf{u} = 0.
```

In two dimensions:

```math
\mathbf{u}(x,y,t) = (u(x,y,t), v(x,y,t)).
```

Thus:

```math
\frac{\partial u}{\partial x}
+
\frac{\partial v}{\partial y}
=
0.
```

### Streamfunction formulation

For two-dimensional incompressible flows, a streamfunction $\psi$ can be introduced so that:

```math
u =
\frac{\partial \psi}{\partial y},
\qquad
v =
-\frac{\partial \psi}{\partial x}.
```

This automatically satisfies incompressibility.

### Vorticity

The scalar vorticity is defined as:

```math
\omega_z =
\frac{\partial v}{\partial x}
-
\frac{\partial u}{\partial y}.
```

Using the streamfunction definition:

```math
\omega_z = -\nabla^2 \psi,
\qquad
\nabla^2 \psi = -\omega_z.
```

### Rate-of-deformation tensor

The symmetric rate-of-deformation tensor is:

```math
\mathbf{D}
=
\frac{1}{2}
\left(
\nabla\mathbf{u}
+
\nabla\mathbf{u}^{T}
\right).
```

For two-dimensional velocity fields:

```math
\mathbf{D}
=
\begin{bmatrix}
\frac{\partial u}{\partial x}
&
\frac{1}{2}
\left(
\frac{\partial u}{\partial y}
+
\frac{\partial v}{\partial x}
\right)
\\[4pt]
\frac{1}{2}
\left(
\frac{\partial u}{\partial y}
+
\frac{\partial v}{\partial x}
\right)
&
\frac{\partial v}{\partial y}
\end{bmatrix}.
```

### Generic viscoelastic constitutive structure

A broad representation of the viscoelastic stress equation may be written as:

```math
\mathrm{Wi}\,\mathcal{D}\boldsymbol{\tau}
+
\mathcal{F}
\left(
\boldsymbol{\tau},
\nabla \mathbf{u}
\right)
=
2(1-\beta)\mathbf{D}.
```

Here:

- $\boldsymbol{\tau}$ is the polymeric extra-stress tensor;
- $\mathrm{Wi}$ is the Weissenberg number;
- $\beta$ is the solvent-viscosity ratio;
- $\mathcal{D}$ denotes a suitable material or objective derivative;
- $\mathcal{F}$ depends on the chosen constitutive model.

The thesis framework involves models such as:

- UCM;
- Oldroyd-B;
- Giesekus;
- LPTT.

### Stress tensor components

```math
\boldsymbol{\tau}
=
\begin{bmatrix}
T_{xx} & T_{xy}
\\
T_{xy} & T_{yy}
\end{bmatrix}.
```

---

## Method of Manufactured Solutions

The **Method of Manufactured Solutions** is a rigorous technique for code verification. Instead of searching for a natural exact solution of the governing equations, one prescribes an artificial but smooth exact solution and computes the forcing terms that make it satisfy the equations exactly.

The manufactured state vector is:

```math
\mathbf{q}_e(x,y,t)
=
\left(
u_e,
v_e,
\psi_e,
\omega_{z,e},
T_{xx,e},
T_{xy,e},
T_{yy,e}
\right).
```

Substituting this manufactured field into the continuous governing equations gives:

```math
\mathcal{L}
\left(
\mathbf{q}_e
\right)
=
\mathbf{f}_{\mathrm{MMS}}.
```

The manufactured source vector contains terms such as:

```math
\mathbf{f}_{\mathrm{MMS}}
=
\left(
f_{\omega_z},
f_{T_{xx}},
f_{T_{xy}},
f_{T_{yy}}
\right).
```

In the original Fortran Case 12 implementation, these are represented by:

```text
tfwz
tfTxx
tfTxy
tfTyy
```

The Python port keeps the same scientific meaning while using NumPy array convention:

```text
(ny, nx)
```

---

## Numerical verification

Given a numerical approximation $\mathbf{q}_h$ on a grid of characteristic spacing $h$, the error is:

```math
\mathbf{e}_h
=
\mathbf{q}_h
-
\mathbf{q}_e.
```

Typical norms include:

```math
\lVert e_h \rVert_{2}
=
\left(
\frac{1}{N}
\sum_{i=1}^{N}
e_i^2
\right)^{1/2},
```

and:

```math
\lVert e_h \rVert_{\infty}
=
\max_i |e_i|.
```

The observed order of convergence is estimated by:

```math
p
=
\frac{
\log(E_{h1}/E_{h2})
}{
\log(h1/h2)
}.
```

For a fourth-order method, the expected behavior is approximately:

```math
E_h = C h^4,
\qquad
p \approx 4.
```

---

## Thesis figures in this README

The figures below come from the thesis/slides material. The original files were mostly EPS, converted to PNG only so GitHub can render them. The original thesis assets remain preserved separately.

Documentation for the selected figures is available in:

```text
docs/thesis-figures-used-in-readme.md
```

### Visual intuition

<p align="center">
  <img src="assets/thesis/case12/thesis_streamlines.png" width="70%" alt="Streamlines from thesis material">
</p>

### MMS validation pipeline

<p align="center">
  <img src="assets/preview/mms_validation_pipeline.svg" width="90%" alt="MMS validation pipeline">
</p>

### Case 12 exact maps

| $\psi$ | $u$ |
|---|---|
| <img src="assets/thesis/case12/exact_map_psi.png" width="420"> | <img src="assets/thesis/case12/exact_map_u.png" width="420"> |

| $v$ | $\omega_z$ |
|---|---|
| <img src="assets/thesis/case12/exact_map_v.png" width="420"> | <img src="assets/thesis/case12/exact_map_wz.png" width="420"> |

| $T_{xx}$ | $T_{xy}$ | $T_{yy}$ |
|---|---|---|
| <img src="assets/thesis/case12/exact_map_txx.png" width="300"> | <img src="assets/thesis/case12/exact_map_txy.png" width="300"> | <img src="assets/thesis/case12/exact_map_tyy.png" width="300"> |

### Case 12 exact surfaces

| $\psi$ | $u$ |
|---|---|
| <img src="assets/thesis/case12/exact_surface_psi.png" width="420"> | <img src="assets/thesis/case12/exact_surface_u.png" width="420"> |

| $\omega_z$ | $T_{xx}$ |
|---|---|
| <img src="assets/thesis/case12/exact_surface_wz.png" width="420"> | <img src="assets/thesis/case12/exact_surface_txx.png" width="420"> |

### Convergence curves from thesis material

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

### Repository architecture and current source-term status

<p align="center">
  <img src="assets/preview/repository_architecture.svg" width="90%" alt="Repository architecture">
</p>

<p align="center">
  <img src="assets/preview/case12_source_terms_status.svg" width="90%" alt="Case 12 source-term status">
</p>

---

## Current Python status

The Python modernization currently includes:

- package foundation under `src/mms_viscoelastic`;
- grid utilities;
- manufactured solution utilities;
- numerical derivative and error utilities;
- visualization utilities;
- ML-ready dataset construction;
- source-term ports for Case 12;
- Fortran-reference validation for several source terms;
- regression tests with `pytest`;
- examples for reproducible execution.

### Implemented source terms

```text
tfwz   ported and validated against Fortran reference
tfTxx  ported and validated against Fortran reference
tfTxy  ported and validated against Fortran reference
tfTyy  ported; Fortran-reference validation pending
```

### Pending work

```text
tfTyy reference validation
complete Case 12 source-term validation report
solver-level Fortran × Python comparisons
expanded datasets for ML workflows
scientific notebooks
presentation material for thesis defense
```

---

## Case 12 status

Case 12 is currently the main Python modernization target.

### Original reference

```text
src/fortran_original/4Order_LPOG_Euler/case_12.f
```

### Python modules

```text
src/mms_viscoelastic/manufactured_solutions/case_12.py
src/mms_viscoelastic/source_terms/case_12.py
src/mms_viscoelastic/validation/source_term_reference.py
```

### Exact fields

```text
u
v
psi
wz
txx
txy
tyy
```

### Source-term API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import evaluate_case_12_source_terms

grid = Grid2D.from_bounds(nx=65, ny=65)
terms = evaluate_case_12_source_terms(grid)

terms.tfwz
terms.tfTxx
terms.tfTxy
terms.tfTyy
```

---

## Repository structure

```text
.
├── src/
│   ├── fortran_original/
│   │   └── 4Order_LPOG_Euler/
│   │       └── case_12.f
│   │
│   ├── mathematica_original/
│   │
│   └── mms_viscoelastic/
│       ├── grids.py
│       ├── fields.py
│       ├── parameters.py
│       ├── manufactured_solutions/
│       ├── source_terms/
│       ├── validation/
│       ├── visualization/
│       ├── ml/
│       └── numerics/
│
├── thesis/
│   └── original/
│
├── papers/
│   └── original/
│
├── validation/
│   ├── fortran_reference/
│   └── reference_data/
│
├── examples/
├── tests/
├── docs/
├── assets/
│   ├── preview/
│   └── thesis/
│       └── case12/
└── outputs/
```

The following folders are scientific reference material and should not be edited casually:

```text
src/fortran_original/
src/mathematica_original/
thesis/original/
papers/original/
```

---

## Installation

Recommended PowerShell workflow:

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

### Compile check

```powershell
python -m compileall src\mms_viscoelastic
```

### Full test suite

```powershell
python -m pytest -q
```

### Source-term reference validation

```powershell
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```

Expected result:

```text
passed:        True
```

Current validation status:

```text
tfwz   validated
tfTxx  validated
tfTxy  validated
tfTyy  pending validation
```

---

## Examples

### Smoke test

```powershell
python examples\01_smoke_python_package.py
```

### Case 12 manufactured solution

```powershell
python examples\02_case12_manufactured_solution.py
```

### Case 12 visualization and ML dataset

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

### Case 12 source terms

```powershell
python examples\04_case12_tfwz_source_term.py
```

This exports:

```text
outputs/case12/source_terms/case12_ported_sources.csv
```

with columns:

```text
x
y
tfwz
tfTxx
tfTxy
tfTyy
```

### Case 12 reference validation

```powershell
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```

---

## Generated outputs

Generated outputs should remain outside version control:

```text
outputs/
```

Curated figures for the README should be copied to:

```text
assets/preview/
assets/thesis/
```

Recommended rule:

```text
outputs/        → generated, not committed
assets/preview/ → selected diagrams, committed
assets/thesis/  → selected thesis figures, committed
```

---

## Machine Learning preparation

The package includes dataset utilities under:

```text
src/mms_viscoelastic/ml/
```

The goal is to prepare pointwise data from manufactured solutions and source terms for future supervised learning, surrogate modeling and Scientific Machine Learning.

A typical dataset row may include:

```text
x
y
time
beta
u
v
psi
wz
txx
txy
tyy
tfwz
tfTxx
tfTxy
tfTyy
```

Potential future ML tasks:

- learn source-term corrections;
- train surrogate models for manufactured fields;
- build error estimators;
- classify numerical regimes;
- compare physics-informed and data-driven approximations;
- generate reduced-order representations.

---

## Local Git workflow

This project is currently maintained with a local branch/merge/tag workflow.

### Start a new phase

```powershell
git checkout main
git status
git checkout -b feat/name-of-phase
```

### Validate

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
```

### Commit

```powershell
git status
git diff --stat
git add -A
git commit -m "feat: objective message"
```

### Merge locally

```powershell
git checkout main
git merge --no-ff feat/name-of-phase -m "merge: objective message"
```

### Tag

```powershell
git tag -a v0.X.Y -m "v0.X.Y: release summary"
```

### Synchronize later

```powershell
git push origin main
git push origin --tags
```

---

## GitHub About

Suggested short description:

```text
Verification and Python modernization of high-order viscoelastic flow solvers using MMS.
```

Suggested topics:

```text
computational-fluid-dynamics
viscoelastic-flows
non-newtonian-fluids
manufactured-solutions
method-of-manufactured-solutions
mms
oldroyd-b
giesekus
lptt
ucm
finite-differences
compact-schemes
fortran
python
scientific-computing
numerical-methods
numerical-verification
machine-learning
scientific-machine-learning
applied-mathematics
doctoral-thesis
usp
```

---

## Citation

If this repository supports your work, cite the associated doctoral thesis and this repository.

Recommended citation text:

```text
Organista, J. Codes and numerical verification material associated with a doctoral thesis on high-order methods for incompressible viscoelastic flows using the Method of Manufactured Solutions.
```

A `CITATION.cff` file should be maintained as the project approaches a stable release.

---

## Author

**Juniormar Organista**

Applied and Computational Mathematics, computational fluid mechanics, viscoelastic flows, high-order finite differences, Method of Manufactured Solutions and scientific software modernization.

---

## README maintenance rule

This README must not be reduced in future phases.

New phases may add sections, update status tables and improve explanations, but should not remove:

- mathematical formulation;
- MMS explanation;
- verification equations;
- thesis figures;
- installation commands;
- validation commands;
- Git workflow;
- roadmap;
- project structure.

If a section becomes too long, move the expanded discussion to `docs/` and keep a summarized version here with a link. Never replace the README with a minimal placeholder.
