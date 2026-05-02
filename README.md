# Thesis USP Codes — MMS Verification and Python Modernization

Repository associated with the doctoral thesis of **Juniormar Organista**, focused on the numerical verification of high-order solvers for incompressible viscoelastic flows using the **Method of Manufactured Solutions (MMS)**.

The project preserves the original thesis material, Fortran codes and LaTeX sources while progressively building a modern Python package for reproducibility, visualization, dataset generation and future machine learning workflows.

---

## Scientific scope

This repository is organized around four scientific pillars:

1. **Numerical verification**
   - Method of Manufactured Solutions;
   - mesh refinement studies;
   - error norms and convergence order estimation.

2. **Viscoelastic flow modeling**
   - incompressible flows;
   - vorticity-streamfunction formulation;
   - non-Newtonian constitutive models;
   - log-conformation-related variables.

3. **High-order numerical methods**
   - compact finite differences;
   - Poisson solvers;
   - time integration procedures;
   - validation of differential operators.

4. **Python modernization**
   - structured Python package;
   - reproducible examples;
   - tests with `pytest`;
   - ML-ready tabular datasets.

---

## Repository structure

```text
.
├── src/
│   ├── fortran_original/          # Original Fortran sources, preserved
│   ├── mathematica_original/      # Original Mathematica notebooks, preserved
│   └── mms_viscoelastic/          # Python package under construction
│
├── thesis/
│   └── original/                  # Original LaTeX thesis files, preserved
│
├── papers/
│   └── original/                  # Original paper-related files, preserved
│
├── docs/                          # Technical documentation
├── examples/                      # Executable Python examples
├── tests/                         # Python test suite
├── assets/preview/                # Figures used in README and GitHub preview
├── scripts/                       # Auxiliary project scripts
└── outputs/                       # Generated results, ignored by Git
```

---

## Python package

The Python package is located at:

```text
src/mms_viscoelastic/
```

Current components include:

```text
grids.py                         # 2D grid utilities
parameters.py                    # physical/numerical parameters
fields.py                        # field containers
numerics/                        # numerical derivatives and errors
manufactured_solutions/          # manufactured exact solutions
visualization/                   # plotting utilities
ml/                              # ML-ready dataset utilities
```

The package is intentionally modular. The long-term goal is to support:

- Fortran-to-Python validation;
- manufactured-solution generation;
- numerical error analysis;
- scientific visualization;
- dataset generation for machine learning;
- future surrogate modeling and Scientific Machine Learning workflows.

---

## Case 12 — preliminary Python port

The current Python modernization starts with **MMS Case 12**.

Implemented:

- exact manufactured fields;
- field visualization;
- tabular dataset generation;
- train/validation split;
- tests for shape consistency and data integrity.

Generated fields:

```text
u, v, psi, wz, txx, txy, tyy
```

ML-ready columns include coordinates, physical/numerical metadata and engineered features:

```text
x, y, time, beta, x2, y2, xy,
sin_pi_x, cos_pi_x, sin_pi_y, cos_pi_y,
u, v, psi, wz, txx, txy, tyy
```

---

## Preview

![Case 12 manufactured fields](assets/preview/case12_panel.png)

---

## Installation

Recommended workflow on Windows PowerShell:

```powershell
cd C:\dev\Thesis_USP_Codes_Organista_J

python -m venv .venv
.\.venv\Scripts\Activate.ps1

python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements-dev.txt
python -m pip install -e .
```

If PowerShell blocks environment activation:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
.\.venv\Scripts\Activate.ps1
```

---

## Running tests

```powershell
python -m pytest -q
```

Compile-check the package:

```powershell
python -m compileall src\mms_viscoelastic
```

---

## Running examples

Basic package smoke test:

```powershell
python examples\01_smoke_python_package.py
```

Evaluate MMS Case 12:

```powershell
python examples\02_case12_manufactured_solution.py
```

Generate Case 12 figures and datasets:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

Use a finer grid:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py --nx 129 --ny 129
```

Generated outputs are written to:

```text
outputs/case12/
```

---

## Fortran sources

The original Fortran implementation is preserved under:

```text
src/fortran_original/
```

These files are treated as historical and scientific reference material. They should not be modified during repository organization or Python modernization.

Python modules will be built progressively from the Fortran implementation, with tests and documentation added at each step.

---

## LaTeX thesis sources

The original thesis files are preserved under:

```text
thesis/original/
```

The project may later include improved presentation material and documentation, but the original thesis files remain unchanged.

---

## Local Git workflow

This project currently follows a local-maintainer workflow:

```text
branch → changes → tests → commit → local merge into main → local tag → later push to GitHub
```

Recommended sequence:

```powershell
git checkout main
git pull --ff-only origin main

git checkout -b feat/name-of-phase

# make changes
python -m pytest -q

git add .
git commit -m "feat: objective commit message"

git checkout main
git merge --no-ff feat/name-of-phase -m "merge: objective merge message"

git tag -a v0.X.0 -m "v0.X.0: release summary"
```

Synchronize later:

```powershell
git push origin main
git push origin --tags
```

---

## Current roadmap

```text
v0.1.0 — repository reorganization and Python package foundation
v0.2.0 — GitHub documentation foundation
v0.3.0 — MMS Case 12 exact manufactured solution
v0.4.0 — Case 12 visualization and ML dataset pipeline
v0.5.0 — visual README, About and project presentation
v0.6.0 — source-term porting for Case 12
v0.7.0 — Fortran vs Python validation
v0.8.0 — ML dataset expansion
v0.9.0 — scientific notebooks
v1.0.0 — first reproducible scientific release
```

---

## Citation

If this repository supports your work, cite the associated doctoral thesis and this repository. A `CITATION.cff` file is included to support citation on GitHub.

---

## Author

**Juniormar Organista**

Doctoral research in Applied and Computational Mathematics, with emphasis on computational fluid mechanics, viscoelastic flows, high-order numerical methods and manufactured-solution verification.
