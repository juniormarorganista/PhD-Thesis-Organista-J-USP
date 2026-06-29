# Fortran-to-Python Porting Plan

## Goal

Build a Python library called `mms_viscoelastic`, based on the original
Fortran codes of the thesis, prepared for:

- numerical verification;
- data generation;
- error analysis;
- Fortran/Python comparison;
- future integration with Machine Learning and Scientific Machine Learning.

## Porting philosophy

The port must not be a blind line-by-line copy. The original Fortran is
procedural and oriented toward direct numerical execution. The Python layer
must be modular, testable and extensible.

Therefore, the translation will be conceptual:

```text
original Fortran -> mathematical interpretation -> clean Python API -> tests -> validation
```

## Planned modules vs. actual structure

The plan below (`manufactured/`, `solvers/`) was the original sketch; the
package that was actually built renamed/flattened a few of these
directories along the way (`manufactured/` → `manufactured_solutions/`,
`solvers/` folded into `numerics/`) and dropped modules that were never
needed (`models/newtonian.py`, `giesekus.py`, `lptt.py`: the Fortran source
never implements separate modules for these, only parameterizes UCM/
Oldroyd-B via `alphaG`/`xi`; `io/writers.py`: the Python port consumes
Fortran-exported CSVs rather than producing its own). The actual structure,
as of the completion of Etapas 13-16:

```text
mms_viscoelastic/
├── parameters.py
├── grids.py
├── fields.py
├── manufactured_solutions/
│   ├── __init__.py
│   ├── toy.py
│   └── case_01.py ... case_14.py   (all 14 MMS cases)
├── source_terms/
│   ├── __init__.py
│   └── case_01.py ... case_14.py
├── numerics/
│   ├── derivatives.py     (2nd-order + 4th-order compact operators)
│   ├── filters.py         (Lele compact filter)
│   ├── errors.py
│   ├── convergence.py
│   ├── poisson.py
│   ├── time_stepping.py       (non-log Euler step)
│   ├── time_stepping_log.py   (log-conformation Euler step)
│   └── time_marching.py
├── models/
│   ├── ucm.py
│   ├── oldroyd_b.py
│   └── log_conformation.py
├── validation/
├── visualization/
├── ml/      (exploratory extension, not thesis core)
└── io/
```

## Phases

### Phase 1 - Foundation (done)

- installable package;
- grid objects;
- field objects;
- basic errors;
- field -> table conversion;
- minimal tests.

### Phase 2 - Documentation and Fortran map (done)

- professional README;
- project map;
- Fortran file map;
- porting plan;
- GitHub instructions.

### Phase 3 - First MMS case (done)

`case_12.f` was the first Python case ported
(`src/mms_viscoelastic/manufactured_solutions/case_12.py`,
`source_terms/case_12.py`), with exact-field functions, source terms,
Fortran-reference validation, convergence verification and time-marching
(Etapas 1-10, 13). It remains the most extensively validated case (the one
used in the thesis and the JNNFM paper).

### Phase 3b - Remaining 13 MMS cases (done, Etapa 14)

`case_01.f`-`case_11.f` ported and validated against an isolated Fortran
reference per case (PRs #35-#45); `case_13.f`/`case_14.f` identified as
non-MMS lid-driven-cavity benchmarks and ported as trivial modules for
structural parity only (PR #46, `docs/case-13-14-port.md`).

### Phase 4 - Derivatives and verification (done)

2nd-order finite-difference operators (Etapas 8-10) and the thesis's actual
4th-order compact (Padé-type) operators plus the Lele compact filter (Etapa
13), both validated against isolated Fortran reference programs.

### Phase 4b - Log-conformation reformulation (done, Etapa 15)

`conform.f` (eigendecomposition, `T_to_Psi`/`Psi_to_T`,
`kernel_conformation`) and the explicit-Euler variant of `drv_Psi` ported
and validated against an isolated Fortran reference. Not discussed in the
thesis text; ported for source-code fidelity. The RK4 variant and full
multi-step marching remain unported.

### Phase 5 - Machine Learning data (exploratory, not thesis core)

A pipeline producing tables with columns such as:

```text
x, y, t, Re, Wi, beta, alphaG, epsilon, u, v, psi, wz, Txx, Txy, Tyy
```

exists under `src/mms_viscoelastic/ml/` (`datasets.py`, `features.py`,
`splits.py`). No surrogate model has been trained; this phase was
explicitly descoped from the thesis core (see `.claude/PLANEJAMENTO.md`,
Etapa 11). The data could in principle feed supervised regression, error
models, surrogates, neural operators or regime classifiers, but doing so is
not part of this project's deliverables.

### Phase 6 - 4th-order time-marching upgrade (done, Etapa 17)

The Euler-step time-marching pipeline now uses the same scheme as the
original Fortran: the 4th-order LPOG compact derivatives (Etapa 13), the
Lele compact filter applied once per step, and a newly-ported 4th-order
compact multigrid V-cycle Poisson solver (`poi_psi4a`,
`numerics/poisson_multigrid.py`) -- not the pentadiagonal
`band5_poi`/`banbk5_poi`, which a repository-wide grep showed is never
called anywhere in the original Fortran source (dead code). All five
isolated Fortran reference programs (`case12_euler_step`,
`case12_log_euler_step`, `case12_time_step`, `case12_multistep`,
`case12_full_march`) were rewritten to match and their CSVs regenerated;
every comparison still agrees at machine precision (`<1e-12`, mostly
`~1e-16`). See `docs/case-12-4th-order-time-marching.md` for the full
account, including the finding that `case12_full_march`'s long-horizon
divergence (Etapa 12) persists with the real, dissipative scheme too -- just
at a different point in time, not eliminated by it. The RK4 branches
(`solver_rk4`, `solver_log_rk4`) remain unported.

## Quality rule

Every port must follow the cycle:

```text
read Fortran -> document the function -> implement Python -> test -> compare -> version
```

No heroic late-night transcriptions without tests. That is like solving
Navier-Stokes on a napkin: romantic, but dangerous.
