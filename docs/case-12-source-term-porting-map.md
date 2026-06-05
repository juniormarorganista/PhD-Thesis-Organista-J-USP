# MMS Case 12 - source-term porting map

This document maps the original Fortran routine

```text
src/fortran_original/4Order_LPOG_Euler/case_12.f
```

to the Python source-term layer.

## Original Fortran routines

Case 12 contains two key routines:

```fortran
subroutine exact_solutions_case_12(t)
subroutine source_term_case_12(t)
```

The exact solution routine fills:

```text
uxe
uye
wze
psie
Txxe
Txye
Tyye
```

The source-term routine fills:

```text
tfwz
tfTxx
tfTxy
tfTyy
```

## Python target

The Python source-term interface is located at:

```text
src/mms_viscoelastic/source_terms/case_12.py
```

Current status:

```text
tfwz    portado e validado contra referência Fortran (v0.6.1 / v0.6.2)
tfTxx   portado e validado contra referência Fortran (v0.6.3 / v0.6.4)
tfTxy   portado e validado contra referência Fortran (v0.6.5 / v0.6.6)
tfTyy   portado e validado contra referência Fortran (v0.6.7 / PR #5)
```

Todos os quatro termos fonte do Case 12 estão implementados e validados.
`evaluate_case_12_source_terms()` retorna um `Case12SourceTerms` completo.

## Porting order

Recommended order:

1. Port `tfwz`.
2. Validate shape and finite values.
3. Compare selected grid points against Fortran output.
4. Port `tfTxx`.
5. Port `tfTxy`.
6. Port `tfTyy`.
7. Add full regression tests against reference CSV files.
8. Add convergence verification workflow.

## Fortran-to-Python naming table

| Fortran | Python target | Meaning |
|---|---|---|
| `tfwz(i,j)` | `tfwz[j, i]` | vorticity source term |
| `tfTxx(i,j)` | `tfTxx[j, i]` | tensor source term Txx |
| `tfTxy(i,j)` | `tfTxy[j, i]` | tensor source term Txy |
| `tfTyy(i,j)` | `tfTyy[j, i]` | tensor source term Tyy |

The Python convention is NumPy-style:

```text
(ny, nx)
```

Fortran loops use:

```fortran
do j = 1, jmax
  do i = 1, imax
```

Therefore the safe Python mapping is:

```python
array[j, i]
```

using zero-based indices.

## Required validation data

To validate this port properly, generate reference files from the original Fortran code containing at least:

```text
x
y
tfwz
tfTxx
tfTxy
tfTyy
```

for a small grid, for example:

```text
nx = 5
ny = 5
```

and for a production-like grid, for example:

```text
nx = 65
ny = 65
```

## Acceptance criteria

The source-term port will be considered acceptable when:

```text
max_abs_error < 1e-10
```

for double-precision reference values at selected test grids, or when a clear
roundoff-level tolerance is justified in the documentation.
