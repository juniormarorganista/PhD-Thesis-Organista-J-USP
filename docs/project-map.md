# Project Map

This document describes the current organization of the repository after the
initial reorganization.

## Organization principle

The structure separates four layers:

1. **preserved original sources**: Fortran, Mathematica, thesis LaTeX and
   paper material;
2. **technical documentation**: maps, plans and usage instructions;
3. **new Python layer**: library under development;
4. **tests and examples**: validation of the new Python layer.

## Main directories

```text
src/fortran_original/
```

Contains the original Fortran code. Must not be changed during the
reorganization. Any improvement must be made in a separate layer or in a
specific maintenance branch.

```text
src/mathematica_original/
```

Contains the original Mathematica notebooks used to generate manufactured
solutions and/or symbolic terms.

```text
src/mms_viscoelastic/
```

Contains the new Python library. As of the completion of Etapas 13-16, it
holds:

- physical parameters and grids;
- all 14 MMS manufactured solutions and source terms
  (`manufactured_solutions/`, `source_terms/`);
- differential operators, both 2nd-order and the thesis's 4th-order compact
  (Padé-type) operators with the Lele filter (`numerics/`);
- constitutive models and the log-conformation reformulation (`models/`);
- a Poisson solver and Euler time-marching for Case 12 (`numerics/`);
- error routines and reference-data validation utilities (`validation/`);
- data export and an exploratory ML feature/dataset pipeline (`ml/`, not
  part of the thesis core scope; see `README.md`).

See `docs/fortran-code-map.md` for the file-by-file Fortran-to-Python
mapping and current status.

```text
Thesis/
```

Contains the thesis and its original presentation. The original LaTeX files
must be preserved.

```text
Papers/
```

Contains the original paper material and associated figures.

```text
docs/
```

Contains technical documentation and development plans.

```text
tests/
```

Contains automated tests for the Python layer.

```text
examples/
```

Contains minimal runnable examples.

## Evolution convention

New changes must be made on thematic branches:

```text
chore/...   reorganization, infrastructure and maintenance
docs/...    documentation
feat/...    new features
fix/...     fixes
tests/...   tests
slides/...  thesis presentation
```

## Preservation rule

Original thesis, Fortran and Mathematica files must be treated as reference
material. The Python port must be validated against these files, not replace
them.
