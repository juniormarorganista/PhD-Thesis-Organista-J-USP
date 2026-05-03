# Thesis USP Codes — MMS Verification and Python Modernization

Repository associated with the doctoral thesis of **Juniormar Organista**, focused on the numerical verification of high-order solvers for incompressible viscoelastic flows using the **Method of Manufactured Solutions (MMS)**.

The repository preserves the original thesis/Fortran/Mathematica material and builds a Python package for reproducibility, visualization, datasets and future machine-learning workflows.

## Current Python status

Implemented:

- repository reorganization;
- Python package foundation;
- MMS Case 12 exact manufactured solution;
- Case 12 visualization;
- ML-ready dataset generation;
- `tfwz` source-term port and Fortran-reference validation;
- `tfTxx` source-term port and Fortran-reference validation;
- `tfTxy` source-term port and Fortran-reference validation.

Pending:

- `tfTyy`;
- complete source-term validation;
- complete solver port.

## Source terms

Ported and validated:

```text
tfwz
tfTxx
tfTxy
```

Pending:

```text
tfTyy
```

## Validation commands

```powershell
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```

Run all tests:

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
```

## Roadmap

```text
v0.6.6 — validate tfTxy against Fortran reference
v0.6.7 — port tfTyy
v0.6.8 — validate tfTyy against Fortran reference
v0.7.0 — complete Case 12 source-term validation
v1.0.0 — first reproducible scientific release
```
