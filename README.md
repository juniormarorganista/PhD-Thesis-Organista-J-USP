# Thesis USP Codes — MMS Verification and Python Modernization

> **Verificação numérica, transcrição Python e preparação para Machine Learning de códigos de alta ordem para escoamentos viscoelásticos usando o Método das Soluções Manufaturadas.**

Repository associated with the doctoral thesis of **Juniormar Organista**, focused on the numerical verification of high-order solvers for incompressible viscoelastic flows using the **Method of Manufactured Solutions (MMS)**.

Este repositório preserva o material original da tese — códigos Fortran, notebooks Mathematica, arquivos LaTeX e materiais associados — e constrói uma camada moderna em Python para:

- reprodutibilidade científica;
- visualização dos campos manufaturados;
- validação Fortran × Python;
- geração de datasets tabulares;
- preparação futura para aprendizado de máquina científico.

---

## Table of contents

1. [Project vision](#project-vision)
2. [Scientific context](#scientific-context)
3. [Mathematical formulation](#mathematical-formulation)
4. [Method of Manufactured Solutions](#method-of-manufactured-solutions)
5. [Numerical verification](#numerical-verification)
6. [Current Python status](#current-python-status)
7. [Case 12 status](#case-12-status)
8. [Figures and visual material](#figures-and-visual-material)
9. [Repository structure](#repository-structure)
10. [Installation](#installation)
11. [Tests and validation](#tests-and-validation)
12. [Examples](#examples)
13. [Generated outputs](#generated-outputs)
14. [Machine Learning preparation](#machine-learning-preparation)
15. [Local Git workflow](#local-git-workflow)
16. [Roadmap](#roadmap)
17. [GitHub About](#github-about)
18. [Citation](#citation)
19. [Author](#author)
20. [README maintenance rule](#readme-maintenance-rule)

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

The Python layer is not a superficial wrapper. It is being organized as a real scientific package with tests, validation data, examples, source-term ports, visualization utilities and future ML-ready datasets.

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

This is a code-verification problem, not merely a simulation problem. The aim is to verify whether the discretized equations and implementation are mathematically consistent with the expected convergence behavior.

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

---

### Streamfunction formulation

For two-dimensional incompressible flows, a streamfunction $\psi$ can be introduced so that:

```math
u =
\frac{\partial \psi}{\partial y},
\qquad
v =
-\frac{\partial \psi}{\partial x}.
```

This automatically satisfies incompressibility:

```math
\frac{\partial}{\partial x}
\left(
\frac{\partial \psi}{\partial y}
\right)
+
\frac{\partial}{\partial y}
\left(
-\frac{\partial \psi}{\partial x}
\right)
=
0.
```

---

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
\omega_z =
-\nabla^2 \psi.
```

Equivalently:

```math
\nabla^2 \psi = -\omega_z.
```

This Poisson-type equation is central in the vorticity-streamfunction formulation.

---

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

---

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

---

### Stress tensor components

In two-dimensional form, the stress tensor is represented by:

```math
\boldsymbol{\tau}
=
\begin{bmatrix}
T_{xx} & T_{xy}
\\
T_{xy} & T_{yy}
\end{bmatrix}.
```

The manufactured fields therefore include:

```math
T_{xx,e}(x,y,t),
\qquad
T_{xy,e}(x,y,t),
\qquad
T_{yy,e}(x,y,t).
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

instead of the Fortran indexing convention.

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
\log(E_h/E_{h/2})
}{
\log(2)
}.
```

For a fourth-order method, the expected behavior is approximately:

```math
E_h
=
C h^4,
```

so that:

```math
p \approx 4.
```

This is the key diagnostic for numerical code verification.

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

The exact manufactured fields include:

```text
u
v
psi
wz
txx
txy
tyy
```

### Source terms

The source-term interface provides:

```python
evaluate_case_12_tfwz(...)
evaluate_case_12_tfTxx(...)
evaluate_case_12_tfTxy(...)
evaluate_case_12_tfTyy(...)
evaluate_case_12_source_terms(...)
```

The full object is:

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

## Figures and visual material

The README uses curated figures under:

```text
assets/preview/
```

Do not commit the entire `outputs/` directory. Generate outputs locally, then copy only selected figures into `assets/preview/`.

### MMS pipeline

![MMS validation pipeline](assets/preview/mms_validation_pipeline.svg)

### Repository architecture

![Repository architecture](assets/preview/repository_architecture.svg)

### Case 12 source-term status

![Case 12 source-term status](assets/preview/case12_source_terms_status.svg)

### Case 12 manufactured fields

Generate the figure:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

Then copy the selected panel:

```powershell
New-Item -ItemType Directory -Force assets\preview
Copy-Item outputs\case12\figures\case12_panel.png assets\preview\case12_panel.png -Force
```

If available, GitHub renders it here:

![Case 12 manufactured fields](assets/preview/case12_panel.png)

### Additional useful figures

After generation, these may also be copied:

```powershell
Copy-Item outputs\case12\figures\case12_u.png assets\preview\case12_u.png -Force
Copy-Item outputs\case12\figures\case12_v.png assets\preview\case12_v.png -Force
Copy-Item outputs\case12\figures\case12_psi.png assets\preview\case12_psi.png -Force
Copy-Item outputs\case12\figures\case12_wz.png assets\preview\case12_wz.png -Force
Copy-Item outputs\case12\figures\case12_txx.png assets\preview\case12_txx.png -Force
Copy-Item outputs\case12\figures\case12_txy.png assets\preview\case12_txy.png -Force
Copy-Item outputs\case12\figures\case12_tyy.png assets\preview\case12_tyy.png -Force
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
│   └── preview/
└── outputs/
```

### Preserved original material

The following folders are scientific reference material and should not be edited casually:

```text
src/fortran_original/
src/mathematica_original/
thesis/original/
papers/original/
```

### Python package

The active Python package is:

```text
src/mms_viscoelastic/
```

This is the only part intended for systematic modernization and extension.

---

## Installation

Recommended PowerShell workflow:

```powershell
cd C:\dev\Thesis_USP_Codes_Organista_J

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

### Current validation status

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
```

Recommended rule:

```text
outputs/        → generated, not committed
assets/preview/ → selected figures, committed
```

This keeps the repository clean while preserving visual documentation.

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

## Roadmap

```text
v0.6.7.3 — restore reference-comparison mean error metrics
v0.6.8   — validate tfTyy against Fortran reference
v0.7.0   — complete Case 12 source-term validation report
v0.7.1   — improve README and thesis visual documentation
v0.8.0   — expanded ML datasets
v0.8.1   — source-term datasets with validated labels
v0.9.0   — scientific notebooks
v0.9.1   — thesis-presentation figures and explanatory diagrams
v1.0.0   — first reproducible scientific release
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
- figure references;
- installation commands;
- validation commands;
- Git workflow;
- roadmap;
- project structure.

If a section becomes too long, move the expanded discussion to `docs/` and keep a summarized version here with a link. Never replace the README with a minimal placeholder.
