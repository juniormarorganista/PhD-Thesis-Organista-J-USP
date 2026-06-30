# Case 12 multigrid Poisson solver Fortran reference

This directory contains a standalone Fortran program that exports reference
values for the FAS (Full Approximation Scheme) multigrid V-cycle Poisson
solver actually used by the original time-marching loop (Etapa 17).

The solver mirrors the original routines:

```text
src/fortran_original/4Order_LPOG_Euler/poisson.f::poi_psi4a
src/fortran_original/4Order_LPOG_Euler/poisson.f::verifica_poi
src/fortran_original/4Order_LPOG_Euler/poisson.f::function_source
src/fortran_original/4Order_LPOG_Euler/poisson.f::sor
src/fortran_original/4Order_LPOG_Euler/poisson.f::defect
src/fortran_original/4Order_LPOG_Euler/poisson.f::restrict
src/fortran_original/4Order_LPOG_Euler/poisson.f::intpol
src/fortran_original/4Order_LPOG_Euler/poisson.f::correction
src/fortran_original/4Order_LPOG_Euler/poisson.f::tridv
```

The original Fortran thesis code is not modified.

## Scope note: this is the real solver, not `band5_poi`/`banbk5_poi`

A repository-wide search (`grep -rn "band5_poi\|banbk5_poi" src/fortran_original/4Order_LPOG_Euler/*.f`)
shows the pentadiagonal solver in `solvers.f` is **never called** anywhere
in the original Fortran source - it is dead code. The Poisson solver
actually invoked by the time-marching loop is `poi_psi4a`
(`call poi_psi4a` in `subroutine loop`, `main.f:182`), a multigrid V-cycle
over a 4th-order compact (Mehrstellen-type) 9-point discretization, smoothed
by line-relaxation SOR (`sor`, one tridiagonal solve per grid line via
`tridv`). This reference program ports `poi_psi4a`, not `band5_poi`.

## Scope

- Generic `nx`, `ny`, `msh` (number of multigrid levels) are accepted as
  command-line arguments; the default (`33 33 4`) matches `par.nn`'s
  `imax = jmax = 33`, `msh = 4` - the only configuration the original
  Fortran code ever actually runs at (`imax`/`jmax`/`msh` are compile-time
  `parameter`s there, not runtime inputs).
- `psi` is used both as the Dirichlet boundary condition (only the
  outermost row/column are read) and as the V-cycle's initial guess for the
  interior (`function_source`'s `fc(1,:,:) = psi`), matching the Fortran
  exactly: boundary rows are pinned by the tridiagonal system in `sor`
  (`a(1)=0, b(1)=1, c(1)=0`) and never reassigned.
- The convergence tolerance (`erromax = 1.0e-6`) and the V-cycle iteration
  cap (`81`) match `poi_psi4a` exactly.

## Generate the input CSV

The input CSV provides `wz` and `psi` from the already-validated Python
Case 12 exact fields:

```powershell
python validation\fortran_reference\case12_poisson_multigrid\generate_input_csv.py `
  validation\reference_data\case12_poisson_multigrid_input_33x33.csv 33 33 0.0
```

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_poisson_multigrid\case12_poisson_multigrid_reference.f90 `
  -o validation\fortran_reference\case12_poisson_multigrid\case12_poisson_multigrid_reference.exe
```

## Generate reference CSV

Default `33 x 33` grid with `4` multigrid levels:

```powershell
.\validation\fortran_reference\case12_poisson_multigrid\case12_poisson_multigrid_reference.exe `
  validation\reference_data\case12_poisson_multigrid_input_33x33.csv `
  validation\reference_data\case12_poisson_multigrid_output_33x33.csv `
  33 33 4
```

## CSV schema

Input (`case12_poisson_multigrid_input_33x33.csv`):

```text
x,y,wz,psi
```

Output (`case12_poisson_multigrid_output_33x33.csv`):

```text
x,y,psi
```

`psi` is the field after the multigrid V-cycle converges (or after 81
V-cycles, whichever comes first); the boundary values equal the input.

## Validate Python against reference

```powershell
python -m pytest tests\test_case12_poisson_multigrid_fortran_reference.py
```
