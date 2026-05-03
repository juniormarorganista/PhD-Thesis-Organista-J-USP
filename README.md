# Thesis USP Codes — MMS Verification and Python Modernization

Repository associated with the doctoral thesis of **Juniormar Organista**, focused on the numerical verification of high-order solvers for incompressible viscoelastic flows using the **Method of Manufactured Solutions (MMS)**.

The repository has two purposes:

1. preserve the original thesis, Fortran and Mathematica material;
2. build a modern Python package for reproducibility, visualization, datasets and future machine-learning workflows.

---

## Scientific problem

The thesis studies the verification of numerical solvers for incompressible viscoelastic flows. These flows are difficult because the velocity, streamfunction, vorticity and polymeric stresses are coupled through nonlinear constitutive equations.

The numerical question is:

> Does the implemented high-order solver reproduce the theoretical order of accuracy when an exact manufactured solution is imposed?

---

## Mathematical structure

```math
\nabla \cdot \mathbf{u} = 0.
```

```math
u = \frac{\partial \psi}{\partial y},
\qquad
v = -\frac{\partial \psi}{\partial x}.
```

```math
\omega_z =
\frac{\partial v}{\partial x}
-
\frac{\partial u}{\partial y},
\qquad
\nabla^2 \psi = -\omega_z.
```

A generic viscoelastic constitutive representation can be written as:

```math
\mathrm{Wi}\,\mathcal{D}\boldsymbol{\tau}
+
\mathcal{F}(\boldsymbol{\tau}, \nabla \mathbf{u})
=
2(1-\beta)\mathbf{D}.
```

The specific form changes according to the model:

- UCM;
- Oldroyd-B;
- Giesekus;
- LPTT.

---

## Method of Manufactured Solutions

The MMS strategy prescribes an exact field:

```math
\mathbf{q}_e(x,y,t)
=
\left(
u_e, v_e, \psi_e, \omega_{z,e},
T_{xx,e}, T_{xy,e}, T_{yy,e}
\right).
```

The exact field is substituted into the governing equations:

```math
\mathcal{L}(\mathbf{q}_e)
=
\mathbf{f}_{\mathrm{MMS}}.
```

The numerical error is:

```math
\mathbf{e}_h
=
\mathbf{q}_h
-
\mathbf{q}_e.
```

The observed order is estimated by:

```math
p
=
\frac{
\log\left(E_h/E_{h/2}\right)
}{
\log(2)
}.
```

---

## Current Python status

Implemented:

- repository reorganization;
- Python package foundation;
- MMS Case 12 exact manufactured solution;
- Case 12 visualization;
- ML-ready dataset generation;
- `tfwz` source-term port for Case 12;
- `tfTxx` source-term port for Case 12;
- tests for package imports, dataset generation and source-term registry.

Pending:

- `tfTxy`;
- `tfTyy`;
- Fortran-vs-Python regression tests for tensor source terms;
- complete solver port.

---

## Source terms

The Case 12 source-term interface is located at:

```text
src/mms_viscoelastic/source_terms/case_12.py
```

Ported:

```text
tfwz
tfTxx
```

Pending:

```text
tfTxy
tfTyy
```

The full evaluator still raises by default:

```python
evaluate_case_12_source_terms()
```

For explicit partial use:

```python
evaluate_case_12_source_terms(grid, allow_partial=True)
```

This returns `tfwz` and `tfTxx`; pending terms remain `None`.

---

## Preview figures

Generate figures:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

If the image below exists, GitHub will render it:

![Case 12 manufactured fields](assets/preview/case12_panel.png)

---

## Installation

```powershell
cd C:\dev\Thesis_USP_Codes_Organista_J

python -m venv .venv
.\.venv\Scripts\Activate.ps1

python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements-dev.txt
python -m pip install -e .
```

---

## Tests

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
```

---

## Examples

```powershell
python examples\01_smoke_python_package.py
python examples\02_case12_manufactured_solution.py
python examples\03_case12_visualization_and_ml_dataset.py
python examples\04_case12_tfwz_source_term.py
```

---

## Roadmap

```text
v0.6.1 — port tfwz
v0.6.2 — validate tfwz against Fortran reference
v0.6.3 — port tfTxx
v0.6.4 — validate tfTxx against Fortran reference
v0.6.5 — port tfTxy
v0.6.6 — port tfTyy
v0.7.0 — Fortran-vs-Python reference validation
v0.8.0 — expanded ML datasets
v0.9.0 — scientific notebooks
v1.0.0 — first reproducible scientific release
```

---

## Local Git workflow

```powershell
git checkout main
git checkout -b feat/name-of-phase

python -m pytest -q

git add -A
git commit -m "feat: objective message"

git checkout main
git merge --no-ff feat/name-of-phase -m "merge: objective message"

git tag -a v0.X.0 -m "v0.X.0: release summary"
```

Synchronize later:

```powershell
git push origin main
git push origin --tags
```

---

## Author

**Juniormar Organista**

Applied and Computational Mathematics, computational fluid mechanics, viscoelastic flows, high-order finite differences and manufactured-solution verification.
