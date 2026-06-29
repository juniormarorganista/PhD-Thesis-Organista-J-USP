# MMS Case 03 - Python port

This document records the port of the original Fortran file:

```text
src/fortran_original/4Order_LPOG_Euler/case_03.f
```

Case 03 is the third of the 13 MMS cases ported after Case 12 (Etapa 14).
Its velocity, vorticity and stream-function fields (`u`, `v`, `wz`, `psi`)
are algebraically identical to Cases 01 and 02 (verified with SymPy: the
difference expands to exactly `0`). The polymeric stress field uses the
same trigonometric ansatz as Case 12
(`(1 - beta) * (a1 + a3 cos(pi x) + a2 sin(pi x)) * (b1 + b3 cos(pi y) +
b2 sin(pi y))`), reusing
`mms_viscoelastic.manufactured_solutions.case_12.Case12Coefficients` since
the field names and `par.nn` defaults (all equal to one) match exactly.

## Ported components

```text
exact_solutions_case_03(t)   -> manufactured_solutions/case_03.py
source_term_case_03(t)       -> source_terms/case_03.py
boundaries_condictions_case_03(t)  -> not ported separately (see Case 01)
```

## Transcription method

`case_03.f`'s source terms are considerably denser than Cases 01/02 (the
`tfTxy` expression alone spans 225 matched parenthesis pairs). Given the
two genuine transcription slips already caught during the Case 02 port
(one in the hand-written Python, one in the hand-written Fortran
reference), Case 03 was ported with no manual retyping of the algebra at
all:

1. A Fortran-to-Fortran mechanical script strips the Fortran77
   continuation prefix and renames intrinsics
   (`dcos`/`Sin`/`dexp` -> `cos`/`sin`/`exp`), preserving every
   parenthesis character-for-character. Used to build the standalone
   reference program
   (`validation/fortran_reference/case03_fields/case03_fields_reference.f90`).
2. A Fortran-to-Python mechanical script performs the equivalent
   transformation into NumPy syntax (`dcos`/`Sin`/`dexp` ->
   `np.cos`/`np.sin`/`np.exp`, `pii` -> `np.pi`), then splits the result
   into top-level additive terms (`term_1`, `term_2`, ...). Those terms
   were copied verbatim into `source_terms/case_03.py`, never retyped.

Both files are therefore independent mechanical transformations of the
same original text; they agree with each other (and with the Python
exact-solution module) to machine precision, which would be extremely
unlikely if either transcription had introduced an error.

## Python API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case03Parameters, evaluate_case_03_exact
from mms_viscoelastic.source_terms import Case03SourceTermParameters, evaluate_case_03_source_terms

grid = Grid2D.from_bounds(nx=33, ny=33)
field = evaluate_case_03_exact(grid, t=0.0, parameters=Case03Parameters())
terms = evaluate_case_03_source_terms(grid, t=0.0, parameters=Case03SourceTermParameters())
```

## Validation status

- All eleven exported quantities validated against an independent
  standalone Fortran reference program
  (`validation/fortran_reference/case03_fields/`), at two grid sizes
  (`5x5`, `9x9`), to `max_abs_error < 1e-10` (observed: machine precision,
  `~1e-14` worst case on `tfwz`).
- Regression tests: `tests/test_case_03_manufactured.py` (shapes, known
  values, identity of `u`/`v`/`wz`/`psi` with Case 01, ML feature table)
  and `tests/test_case03_fields_fortran_reference.py` (Fortran cross-check).
- Example: `python examples/18_case03_fields_and_source_terms.py`.

## Reproducing

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case03_fields\case03_fields_reference.f90 `
  -o validation\fortran_reference\case03_fields\case03_fields_reference.exe
.\validation\fortran_reference\case03_fields\case03_fields_reference.exe `
  validation\reference_data\case03_fields_9x9_default.csv 9 9

python examples\18_case03_fields_and_source_terms.py
python -m pytest tests\test_case_03_manufactured.py tests\test_case03_fields_fortran_reference.py -v
```
