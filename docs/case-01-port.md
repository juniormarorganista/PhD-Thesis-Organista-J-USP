# MMS Case 01 - Python port

This document records the port of the original Fortran file:

```text
src/fortran_original/4Order_LPOG_Euler/case_01.f
```

Case 01 is the first of the 13 MMS cases ported after Case 12 (Etapa 14 in
the local planning file). Unlike Case 12, its manufactured polymeric stress
field (`txx`, `txy`, `tyy`) is spatially uniform: it reproduces the fixed
target values `c1=1`, `c2=2`, `c3=3` (compile-time `parameter`s in
`par.nn`), scaled by `(1 - beta)`. Only the velocity/vorticity/stream
function fields (`u`, `v`, `wz`, `psi`) vary with position, via the
polynomial `(x-1)^2 x^2 (\ldots)` shape used throughout this case.

## Ported components

```text
exact_solutions_case_01(t)   -> manufactured_solutions/case_01.py
source_term_case_01(t)       -> source_terms/case_01.py
boundaries_condictions_case_01(t)  -> not ported separately (see below)
```

`boundaries_condictions_case_01` re-evaluates the same exact-solution
formulas restricted to the boundary nodes (`i = 1, imax` and
`j = 1, jmax`); it has no distinct algebraic expression. The Python port
therefore does not include a separate boundary-condition function for
Case 01, the same convention already used for Case 12 (see
`mms_viscoelastic.numerics.time_marching`, which re-imposes the exact
solution at the boundary rather than reading a dedicated boundary formula).

## Verification of the exact fields

The stream-function identities

```text
u  = d(psi)/dy
v  = -d(psi)/dx
wz = d(u)/dy - d(v)/dx
```

were checked symbolically with SymPy against the transcribed expressions
for `u`, `v`, `psi` and `wz`; all three residuals simplify to exactly `0`,
confirming the transcription is self-consistent (same convention as Case 12:
`wz` equals `+laplacian(psi)`, not its negative).

## Python API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case01Parameters, evaluate_case_01_exact
from mms_viscoelastic.source_terms import Case01SourceTermParameters, evaluate_case_01_source_terms

grid = Grid2D.from_bounds(nx=33, ny=33)
field = evaluate_case_01_exact(grid, t=0.0, parameters=Case01Parameters())
terms = evaluate_case_01_source_terms(grid, t=0.0, parameters=Case01SourceTermParameters())
```

`field.txx`/`field.txy`/`field.tyy` are returned as full `(ny, nx)` arrays
(broadcast from the underlying scalar expression) so that all `Field2D`
components share the same shape, as required by
`Field2D.as_feature_matrix`.

## Validation status

- Symbolic self-consistency of the exact fields (SymPy), see above.
- All eleven exported quantities (`u`, `v`, `wz`, `psi`, `txx`, `txy`,
  `tyy`, `tfwz`, `tfTxx`, `tfTxy`, `tfTyy`) validated against an
  independent standalone Fortran reference program
  (`validation/fortran_reference/case01_fields/`), at two grid sizes
  (`5x5`, `9x9`), to `max_abs_error < 1e-10` (observed: machine precision,
  `~1e-15` to `~1e-17`).
- Regression tests: `tests/test_case_01_manufactured.py` (shapes, known
  values, ML feature table) and
  `tests/test_case01_fields_fortran_reference.py` (Fortran cross-check).
- Example: `python examples/16_case01_fields_and_source_terms.py`.

## Reproducing

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case01_fields\case01_fields_reference.f90 `
  -o validation\fortran_reference\case01_fields\case01_fields_reference.exe
.\validation\fortran_reference\case01_fields\case01_fields_reference.exe `
  validation\reference_data\case01_fields_9x9_default.csv 9 9

python examples\16_case01_fields_and_source_terms.py
python -m pytest tests\test_case_01_manufactured.py tests\test_case01_fields_fortran_reference.py -v
```
