# MMS Cases 13 and 14 - Python port

This document records the port of the original Fortran files:

```text
src/fortran_original/4Order_LPOG_Euler/case_13.f
src/fortran_original/4Order_LPOG_Euler/case_14.f
```

Cases 13 and 14 are the twelfth and thirteenth (and last) of the 13 MMS
cases tracked in Etapa 14 -- but unlike Cases 01-12, **they are not
method-of-manufactured-solutions verification cases**.

## Why these two are different

Reading both Fortran files shows that `exact_solutions_case_13/14(t)`
sets every field (`uxe`, `uye`, `wze`, `psie`, `Txxe`, `Txye`, `Tyye`) to
`0.0` identically, and `source_term_case_13/14(t)` sets every source
term (`tfwz`, `tfTxx`, `tfTxy`, `tfTyy`) to `0.0` identically. There is
no manufactured exact solution and no forcing to derive symbolically --
which is exactly why neither case has a corresponding Mathematica
notebook under `src/mathematica_original/` (every one of Cases 01-12
does).

Instead, both cases are physical benchmark configurations of the
classic **lid-driven cavity flow** problem, driven entirely through
`boundaries_condictions_case_13/14(t)`:

- **Case 13**: a constant unit lid velocity, `ux(i, jmax) = 1.0`.
- **Case 14**: a smooth-startup lid velocity that ramps up in time via a
  hyperbolic tangent, `ux(i, jmax) = 0.8 * (1 - x**2)**2 * x**2 * (1 +
  tanh(8*t - 4))`.

These are used to exercise the solver's qualitative behavior (vortex
formation, startup transients) rather than to verify numerical order of
accuracy -- there is no exact solution to compare against.

## Scope decision

Given the above, porting these two as conventional MMS cases (with a
Fortran cross-validation, term-by-term algebra, etc.) would have no
verification value: both sides of any comparison would be `0.0`. After
confirming this finding, the chosen scope was to port them as **trivial
modules**, matching the same function signatures as every other case
for structural parity, but documenting plainly that they carry no
numerical content. Their only non-trivial content -- the boundary lid
velocity profile -- is **not** ported, following the same convention
already used for every other case (boundary conditions are never ported
separately).

## Ported components

```text
exact_solutions_case_13(t)   -> manufactured_solutions/case_13.py (returns all-zero fields)
source_term_case_13(t)       -> source_terms/case_13.py (returns all-zero source terms)
exact_solutions_case_14(t)   -> manufactured_solutions/case_14.py (returns all-zero fields)
source_term_case_14(t)       -> source_terms/case_14.py (returns all-zero source terms)
boundaries_condictions_case_13/14(t)  -> not ported (same convention as every other case)
```

## Python API

```python
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_case_13_exact, evaluate_case_14_exact
from mms_viscoelastic.source_terms import evaluate_case_13_source_terms, evaluate_case_14_source_terms

grid = Grid2D.from_bounds(nx=33, ny=33)
field_13 = evaluate_case_13_exact(grid, t=0.0)  # every component is 0.0
terms_13 = evaluate_case_13_source_terms(grid, t=0.0)  # every component is 0.0
```

## Validation status

- No Fortran cross-validation was performed (both sides would trivially
  be `0.0`; this would not exercise any transcription risk).
- Regression tests: `tests/test_case_13_manufactured.py` and
  `tests/test_case_14_manufactured.py` (shapes, confirming all fields and
  source terms are identically zero, ML feature table).
- Example: `python examples/27_case13_case14_trivial_benchmarks.py`.

## Reproducing

```powershell
python examples\27_case13_case14_trivial_benchmarks.py
python -m pytest tests\test_case_13_manufactured.py tests\test_case_14_manufactured.py -v
```
