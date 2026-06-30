# Fortran Code Map

Reorganized original directory:

```text
src/fortran_original/4Order_LPOG_Euler/
```

This map tracks the port to Python and its current status. See
`.claude/PLANEJAMENTO.md` (local, not in the repository) for the staged
plan; Etapas 13-17 are complete (Etapa 17 closes the gap: see
`docs/case-12-4th-order-time-marching.md`).

## Main files

| File | Role in the original code | Python destination | Status |
|---|---|---|---|
| `main.f` | main program and simulation orchestration | -- | not ported; the project keeps a library of validated pieces rather than a single runnable equivalent (see Etapa 13 decisions) |
| `comm.var`, `comm.vare` | common variables/global parameters | `mms_viscoelastic.parameters` | ported (`FlowParameters`) |
| `par.nn`, `par.nn.template` | run parameters | `mms_viscoelastic.parameters`, `mms_viscoelastic.numerics.filters` (`alphaf`/`af`/`bf`/`cf`/`df`) | ported |
| `case_12.f` | MMS case 12 (used in the thesis and the JNNFM paper) | `mms_viscoelastic.manufactured_solutions.case_12`, `mms_viscoelastic.source_terms.case_12` | ported and validated against Fortran reference (most extensive validation: source terms, convergence, time-marching, log-conformation) |
| `case_01.f`--`case_11.f` | remaining genuine MMS cases | `mms_viscoelastic.manufactured_solutions.case_NN`, `mms_viscoelastic.source_terms.case_NN` | ported and validated against an isolated Fortran reference per case (Etapa 14, PRs #35-#45) |
| `case_13.f`, `case_14.f` | lid-driven-cavity benchmark configurations (not MMS: exact solution and source term are identically zero, no Mathematica notebook) | `mms_viscoelastic.manufactured_solutions.case_13`/`case_14`, `mms_viscoelastic.source_terms.case_13`/`case_14` | ported as trivial (zeroed) modules for structural parity only, no Fortran cross-validation (both sides trivially zero); see `docs/case-13-14-port.md` (Etapa 14, PR #46) |
| `derivs.f` | compact 4th-order derivatives | `mms_viscoelastic.numerics.derivatives` (`compact_derivative_x/y`, `compact_second_derivative_x/y`) | ported and validated against Fortran reference (Etapa 13); now used throughout the time-marching pipeline (Etapa 17, see below) |
| `drvs.f` | nonlinear/convective terms (`drv`, `nlterms`) and their log-conformation counterparts (`drv_Psi`, `nlterms_logPsi`) | `mms_viscoelastic.numerics.time_stepping` (non-log half), `mms_viscoelastic.numerics.time_stepping_log` (log-conformation half) | both halves ported and validated against Fortran reference using the real 4th-order compact operators (Etapa 17, was a 2nd-order placeholder before); the log-conformation half covers the explicit-Euler variant only (Etapa 15) |
| `poisson.f` | Poisson equation solver (`poi_psi4a`: 4th-order compact 9-point/Mehrstellen stencil + FAS multigrid V-cycle) | `mms_viscoelastic.numerics.poisson_multigrid` (`solve_poisson_multigrid_4th_order`) | ported and validated to exact bit-for-bit agreement against an isolated Fortran reference, plus 4th-order convergence (Etapa 17). A separate 2nd-order 5-point direct solve (`mms_viscoelastic.numerics.poisson`, Etapa 8) remains as an independent utility, no longer used by the time-marching pipeline |
| `solvers.f` | Thomas algorithm (`tridx`/`tridy`) and pentadiagonal Poisson solve (`band5_poi`/`banbk5_poi`) | `mms_viscoelastic.numerics.derivatives` (`solve_tridiagonal`) | `tridx`/`tridy` ported; `band5_poi`/`banbk5_poi` **not ported -- confirmed dead code**: a repository-wide search shows these two subroutines are never called anywhere in the original Fortran source (Etapa 17 finding, see `docs/case-12-4th-order-time-marching.md` Section 1). `poi_psi4a` (above) is the real, actually-called Poisson solver |
| `euler_rk.f` | Euler/RK4 time integration, including the `ty_sim==22` log-conformation Euler branch (`solver_log_euler`) | `mms_viscoelastic.numerics.time_stepping`, `mms_viscoelastic.numerics.time_marching`, `mms_viscoelastic.numerics.time_stepping_log` | Euler branches ported with the real 4th-order compact + Lele filter + multigrid Poisson pipeline (Etapa 17) and validated against Fortran reference; the RK4 branches (`solver_rk4`, `solver_log_rk4`) are not ported |
| `conform.f` | log-conformation transformation (`decomp_matrix`, `multi_matrix`, `T_to_Psi`, `Psi_to_T`, `tf_conformation`, `kernel_conformation`) | `mms_viscoelastic.models.log_conformation` | ported and validated against Fortran reference (Etapa 15); not used in the thesis text - ported for source-code fidelity only. `cc_Psi` (boundary-specific re-projection) not ported |
| `condictions.f` | initial/boundary conditions, MMS case dispatch | `mms_viscoelastic.numerics.time_marching` (boundary re-imposition for case 12 only) | partially ported; the generic `SELECT CASE (ty_mms)` dispatch across all 14 cases has no Python equivalent |
| `boundarie_vorticity_cavity_flow.f` | vorticity boundary conditions for the cavity problem | -- | not ported; case 12's time-marching re-imposes the exact manufactured boundary condition instead (a documented deviation, see `docs/case-12-time-marching-comparison.md`) |
| `normerror.f` | error norm computation | `mms_viscoelastic.numerics.errors` | ported |
| `write_solutions.f` | results output | `mms_viscoelastic.io` | ported (`read_reference_csv`; writing is not needed since the Python port consumes Fortran-exported CSVs rather than producing its own) |
| `loop.f` | main simulation loops, Lele compact filter | `mms_viscoelastic.numerics.filters` (`lele_filter_x/y`) | filter ported and validated against Fortran reference (Etapa 13); now applied once per time step by `mms_viscoelastic.numerics.time_marching` (Etapa 17), matching `subroutine loop`'s `filter_trid_x`/`filter_trid_y`/`poi_psi4a`/`calc_ux`/`calc_uy` order. The simulation-loop orchestration itself (the `do` loop calling `loop` every step in `main.f`) has no single Python equivalent, only the per-step pieces |
| `formatted.f` | auxiliary formatting/output | -- | not ported; not needed (see `write_solutions.f`) |
| `makefile` | Fortran build | keep in `src/fortran_original` | n/a |

## Porting order followed so far

1. `comm.var`, `comm.vare`, `par.nn` -> Python parameters (done);
2. `case_12.f` -> first reference MMS case (done);
3. `normerror.f` -> errors and convergence order (done);
4. `euler_rk.f`, `solvers.f` (`tridx`/`tridy`), `poisson.f` (2nd-order) ->
   time-marching pipeline (done, 2nd-order placeholder, since superseded);
5. `derivs.f`, `loop.f` (filter) -> compact 4th-order derivatives and Lele
   filter (done, Etapa 13);
6. remaining cases `case_01.f`--`case_11.f`, `case_13.f`, `case_14.f`
   (done, Etapa 14, PRs #35-#46);
7. `conform.f`, log-conformation half of `drvs.f` (done, Etapa 15, Euler
   variant only; RK4 variant deferred);
8. time-marching pipeline upgrade to the real scheme: `poisson.f`'s
   `poi_psi4a` (4th-order compact multigrid V-cycle -- not the pentadiagonal
   `band5_poi`/`banbk5_poi`, confirmed dead code), and reapointing
   `time_stepping.py`/`time_stepping_log.py`/`time_marching.py` to the
   compact derivatives and Lele filter (done, Etapa 17, see
   `docs/case-12-4th-order-time-marching.md`).

## Module acceptance criteria

A ported module should only be considered accepted when it:

- has a unit test;
- has a minimal example;
- has short documentation;
- reproduces a known output or an expected mathematical property;
- does not change the original Fortran files.
