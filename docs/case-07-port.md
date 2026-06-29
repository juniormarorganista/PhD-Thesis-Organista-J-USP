# MMS Case 07 - Python port

This document records the port of the original Fortran file:

```text
src/fortran_original/4Order_LPOG_Euler/case_07.f
```

Case 07 is the seventh of the 13 MMS cases ported after Case 12 (Etapa
14). Unlike every other ported case, its velocity, vorticity and
stream-function fields (`u`, `v`, `wz`, `psi`) decay in time through a
hyperbolic tangent, `-1 + tanh(4 - 8*at*t)`, rather than `exp(-at*t)`.
The polymeric stress field is spatially uniform, algebraically identical
to Case 01's ansatz.

The `u`/`v`/`wz`/`psi` formulas were verified self-consistent via the
stream-function identities (`u = d(psi)/dy`, `v = -d(psi)/dx`,
`wz = Laplacian(psi)`, matching the sign convention already established
for Case 01) by hand, with SymPy-equivalent algebra confirming exact
agreement before the Fortran cross-validation ran.

## Ported components

```text
exact_solutions_case_07(t)   -> manufactured_solutions/case_07.py
source_term_case_07(t)       -> source_terms/case_07.py
boundaries_condictions_case_07(t)  -> not ported separately (see Case 01)
```

## New transcription cases

`case_07.f` introduced two patterns the mechanical scripts had not
handled before:

- Hyperbolic functions `Tanh`/`Cosh` (the time-derivative of `tanh` is
  `1/cosh**2`, appearing in `tfwz`). Added to both mechanical scripts'
  intrinsic-rename tables (`np.tanh`/`np.cosh` for the Python script;
  Fortran's generic `tanh`/`cosh` already work without renaming).
- Fortran double-precision literal exponents (`8.0d0`, `1.d0`). The
  Fortran-to-Python script now converts the `d`/`D` exponent marker to
  `e` whenever it immediately follows a digit or decimal point, so
  `8.0d0` becomes the valid Python literal `8.0e0` without disturbing
  `dexp`/`dcos`/etc. (whose leading `d` starts at a word boundary, never
  preceded by a digit).

## A floating-point evaluation-order artifact in the Fortran reference

Cross-validating against the standalone Fortran reference initially
showed `tfTyy` disagreeing by ~4.5e-8 -- far above the `<1e-10` target
and far above ordinary double-precision rounding (~1e-16), while the
other ten quantities matched at machine precision immediately. Given
the close textual match already established for every other case, this
warranted direct investigation rather than assuming a transcription slip:

1. A byte-for-byte comparison confirmed the compiled Fortran statement
   for `tfTyy` is identical (modulo whitespace) to the text mechanically
   extracted from `case_07.f`, and that text is in turn algebraically
   identical to the Python module's `evaluate_case_07_tftyy`.
2. Evaluating the same text in Python/NumPy reproduced the
   mathematically expected value (~1.5 at the grid origin, exact to
   ~1e-16), matching a Fortran re-implementation that assigns the same
   sub-expressions to named intermediate variables before combining
   them.
3. The discrepancy reproduced under both `-O0` and `-O2`, and under
   `-ffp-contract=off`, ruling out compiler optimization level and FMA
   contraction as the cause.

The conclusion: gfortran's evaluation of this one deeply nested,
single-statement expression exhibits a floating-point cancellation
artifact (the numerator and denominator both involve `exp(8)`-scale
terms that nearly cancel to an O(1) result), not present when the same
algebra is split into named sub-expressions first. The Fortran reference
program's `tftyy_case07` therefore uses the decomposed form -- the same
style already used throughout the Python port's `term_1`, `term_2`,
... decomposition, applied here for numerical robustness rather than
transcription safety.

## Python API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case07Parameters, evaluate_case_07_exact
from mms_viscoelastic.source_terms import Case07SourceTermParameters, evaluate_case_07_source_terms

grid = Grid2D.from_bounds(nx=33, ny=33)
field = evaluate_case_07_exact(grid, t=0.0, parameters=Case07Parameters())
terms = evaluate_case_07_source_terms(grid, t=0.0, parameters=Case07SourceTermParameters())
```

## Validation status

- All eleven exported quantities validated against an independent
  standalone Fortran reference program
  (`validation/fortran_reference/case07_fields/`), at two grid sizes
  (`5x5`, `9x9`), to `max_abs_error < 1e-10` (observed: machine precision
  after decomposing `tftyy_case07`).
- Regression tests: `tests/test_case_07_manufactured.py` (shapes, known
  values, identity of `txx`/`txy`/`tyy` with Case 01, ML feature table)
  and `tests/test_case07_fields_fortran_reference.py` (Fortran cross-check).
- Example: `python examples/22_case07_fields_and_source_terms.py`.

## Reproducing

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case07_fields\case07_fields_reference.f90 `
  -o validation\fortran_reference\case07_fields\case07_fields_reference.exe
.\validation\fortran_reference\case07_fields\case07_fields_reference.exe `
  validation\reference_data\case07_fields_9x9_default.csv 9 9

python examples\22_case07_fields_and_source_terms.py
python -m pytest tests\test_case_07_manufactured.py tests\test_case07_fields_fortran_reference.py -v
```
