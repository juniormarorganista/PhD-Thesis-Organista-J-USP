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

Contains the new Python library. This directory will be expanded to hold:

- physical parameters;
- grids;
- numerical fields;
- manufactured solutions;
- differential operators;
- solvers;
- error routines;
- data export;
- Machine Learning integration.

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
