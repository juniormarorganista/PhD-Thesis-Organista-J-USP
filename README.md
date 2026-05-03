# Thesis USP Codes — MMS Verification and Python Modernization

Repository associated with the doctoral thesis of **Juniormar Organista**, focused on numerical verification of high-order solvers for incompressible viscoelastic flows using the **Method of Manufactured Solutions (MMS)**.

## Current Python status

Implemented:

- MMS Case 12 exact manufactured solution;
- Case 12 visualization;
- ML-ready dataset generation;
- `tfwz`, `tfTxx`, `tfTxy` source-term ports and Fortran-reference validation;
- `tfTyy` source-term port;
- validation compatibility with the current `Grid2D.from_bounds(nx, ny)` API;
- reference-comparison metadata fields `nx`, `ny` and `n_points`.

Pending:

- Fortran-reference validation for `tfTyy`;
- complete source-term validation report;
- complete solver port.

## Main validation commands

```powershell
python -m compileall src\mms_viscoelastic
python -m pytest -q
python examples\05_validate_case12_tfwz_reference.py
python examples\06_validate_case12_tfTxx_reference.py
python examples\07_validate_case12_tfTxy_reference.py
```

## Roadmap

```text
v0.6.7.2 — restore reference-comparison grid metadata
v0.6.8   — validate tfTyy against Fortran reference
v0.7.0   — complete Case 12 source-term validation
v1.0.0   — first reproducible scientific release
```
