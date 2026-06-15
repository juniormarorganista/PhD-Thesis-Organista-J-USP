# Initial Repository Map

This repository brings together the original thesis files, Fortran codes,
Mathematica notebooks, the scientific paper and the new Python foundation.

## Preservation rule

The original Fortran codes and LaTeX files must not be changed at this stage.
They have only been reorganized with `git mv`, preserving the Git history.

## Structure after reorganization

```text
src/
  fortran_original/        # original thesis Fortran code
  mathematica_original/    # original Mathematica notebooks
  mms_viscoelastic/        # new Python library under construction

thesis/
  original/                # original LaTeX thesis

papers/
  original/                # original paper and figures

docs/                      # repository technical documentation
examples/                  # runnable examples for the Python library
tests/                     # automated tests for the Python library
```

## Goal of the Python foundation

The `src/mms_viscoelastic` folder was created to evolve into a Python library.
The Fortran -> Python port must be incremental, testable and driven by
comparison against the original implementation.
