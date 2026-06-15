# Initial Fortran Code Map

Reorganized original directory:

```text
src/fortran_original/4Order_LPOG_Euler/
```

This map guides the future port to Python.

## Main files

| File | Role in the original code | Planned Python destination |
|---|---|---|
| `main.f` | main program and simulation orchestration | `mms_viscoelastic.runners` |
| `comm.var`, `comm.vare` | common variables/global parameters | `mms_viscoelastic.parameters` |
| `par.nn`, `par.nn.template` | run parameters | `mms_viscoelastic.config` |
| `case_01.f`--`case_14.f` | MMS cases and source terms | `mms_viscoelastic.manufactured.cases` |
| `derivs.f` | derivatives and differential operators | `mms_viscoelastic.numerics.derivatives` |
| `drvs.f` | possible derivative variant/helper | evaluate before porting |
| `poisson.f` | Poisson equation solver | `mms_viscoelastic.solvers.poisson` |
| `solvers.f` | auxiliary solver routines | `mms_viscoelastic.solvers.linear` |
| `euler_rk.f` | Euler/RK time integration | `mms_viscoelastic.solvers.time` |
| `conform.f` | log-conformation transformation | `mms_viscoelastic.models.log_conformation` |
| `condictions.f` | initial/boundary conditions | `mms_viscoelastic.boundary` |
| `boundarie_vorticity_cavity_flow.f` | vorticity boundary conditions for the cavity problem | `mms_viscoelastic.boundary.vorticity` |
| `normerror.f` | error norm computation | `mms_viscoelastic.numerics.errors` |
| `write_solutions.f` | results output | `mms_viscoelastic.io` |
| `loop.f` | main simulation loops | `mms_viscoelastic.runners` |
| `formatted.f` | auxiliary formatting/output | `mms_viscoelastic.io` |
| `makefile` | Fortran build | keep in `src/fortran_original` |

## Recommended porting order

1. `comm.var`, `comm.vare`, `par.nn` -> Python parameters;
2. `case_12.f` -> first reference MMS case;
3. `normerror.f` -> errors and convergence order;
4. `derivs.f` -> second-order derivatives first;
5. `poisson.f` -> simplified Poisson solver;
6. `euler_rk.f` -> time integration;
7. `conform.f` -> log-conformation;
8. remaining cases `case_01.f`--`case_14.f`.

## Module acceptance criteria

A ported module should only be considered accepted when it:

- has a unit test;
- has a minimal example;
- has short documentation;
- reproduces a known output or an expected mathematical property;
- does not change the original Fortran files.
