# MMS Case 04 - Python port

This document records the port of the original Fortran file:

```text
src/fortran_original/4Order_LPOG_Euler/case_04.f
```

Case 04 is the fourth of the 13 MMS cases ported after Case 12 (Etapa
14). Its velocity, vorticity and stream-function fields (`u`, `v`, `wz`,
`psi`) follow a Taylor-Green-vortex trigonometric ansatz
(`u = 8 cos(pi x) sin(pi y) / decay`, etc.), distinct from cases
01/02/03 (polynomial) and case 12 (polynomial-times-trigonometric). The
polymeric stress field is spatially uniform, algebraically identical to
case 01's ansatz (verified directly: same `-((-1+beta)*c_i)/decay`
formula for `i = 1, 2, 3`).

## Ported components

```text
exact_solutions_case_04(t)   -> manufactured_solutions/case_04.py
source_term_case_04(t)       -> source_terms/case_04.py
boundaries_condictions_case_04(t)  -> not ported separately (see Case 01)
```

## A sign error caught by cross-validation

Case 04's source terms are far more compact than Cases 02/03 -- the
largest, `tfTxx`/`tfTyy`, span only 19 matched parenthesis pairs.
Because they read at a glance, they were transcribed by hand into
`source_terms/case_04.py` instead of being copied verbatim from a
mechanically-generated Python expression (the discipline used for Cases
02 and 03).

This reintroduced exactly the kind of error the mechanical-copy
discipline exists to prevent: the outer factor multiplying
`tfTxx`/`tfTxy`/`tfTyy` in the original Fortran is `(-1 + betann)`, but
the hand-written Python used `(1.0 - betann)` -- the negative of the
correct value (`-1 + beta` and `1 - beta` are not the same expression).

Cross-validating against the independent Fortran reference caught it
immediately: `tfTxx`, `tfTxy` and `tfTyy` showed constant absolute
errors of roughly 100, 52 and 198 respectively, the same at every grid
point regardless of `x` or `y` -- the signature of a sign error in a
term that does not depend on position (the `(-1+betann)` factor multiplies
the *entire* expression, including its position-independent pieces).
`u`, `v`, `wz`, `psi`, `txx`, `txy`, `tyy` and `tfwz` already matched to
machine precision, narrowing the search immediately to the three
remaining tensor source terms. Fixed by restoring `(-1.0 + betann)` in
all three functions.

## Python API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case04Parameters, evaluate_case_04_exact
from mms_viscoelastic.source_terms import Case04SourceTermParameters, evaluate_case_04_source_terms

grid = Grid2D.from_bounds(nx=33, ny=33)
field = evaluate_case_04_exact(grid, t=0.0, parameters=Case04Parameters())
terms = evaluate_case_04_source_terms(grid, t=0.0, parameters=Case04SourceTermParameters())
```

## Validation status

- All eleven exported quantities validated against an independent
  standalone Fortran reference program
  (`validation/fortran_reference/case04_fields/`), at two grid sizes
  (`5x5`, `9x9`), to `max_abs_error < 1e-10` (observed: machine
  precision, after fixing the sign error above).
- Regression tests: `tests/test_case_04_manufactured.py` (shapes, known
  values, identity of `txx`/`txy`/`tyy` with Case 01, ML feature table)
  and `tests/test_case04_fields_fortran_reference.py` (Fortran cross-check).
- Example: `python examples/19_case04_fields_and_source_terms.py`.

## Reproducing

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case04_fields\case04_fields_reference.f90 `
  -o validation\fortran_reference\case04_fields\case04_fields_reference.exe
.\validation\fortran_reference\case04_fields\case04_fields_reference.exe `
  validation\reference_data\case04_fields_9x9_default.csv 9 9

python examples\19_case04_fields_and_source_terms.py
python -m pytest tests\test_case_04_manufactured.py tests\test_case04_fields_fortran_reference.py -v
```
