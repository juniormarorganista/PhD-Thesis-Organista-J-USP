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

## Planned modules

```text
mms_viscoelastic/
├── parameters.py
├── grids.py
├── fields.py
├── manufactured/
│   ├── __init__.py
│   ├── base.py
│   ├── case_12.py
│   └── cases.py
├── numerics/
│   ├── derivatives.py
│   ├── errors.py
│   ├── convergence.py
│   └── filters.py
├── models/
│   ├── newtonian.py
│   ├── ucm.py
│   ├── oldroyd_b.py
│   ├── giesekus.py
│   ├── lptt.py
│   └── log_conformation.py
├── solvers/
│   ├── poisson.py
│   ├── linear.py
│   └── time.py
├── ml/
│   ├── datasets.py
│   ├── features.py
│   └── splits.py
└── io/
    ├── readers.py
    └── writers.py
```

## Phases

### Phase 1 - Foundation

Already started:

- installable package;
- grid objects;
- field objects;
- basic errors;
- field -> table conversion;
- minimal tests.

### Phase 2 - Documentation and Fortran map

Current:

- professional README;
- project map;
- Fortran file map;
- porting plan;
- GitHub instructions.

### Phase 3 - First MMS case

Implement `case_12.f` as the first Python case:

```text
src/mms_viscoelastic/manufactured/case_12.py
```

Must contain functions for:

- `u_exact(x, y, t, params)`;
- `v_exact(x, y, t, params)`;
- `psi_exact(x, y, t, params)`;
- `wz_exact(x, y, t, params)`;
- `stress_exact(x, y, t, params)`;
- source terms once ported.

### Phase 4 - Derivatives and verification

Implement finite-difference derivatives and test them using smooth analytic
functions.

### Phase 5 - Machine Learning data

Create a pipeline that produces tables with columns such as:

```text
x, y, t, Re, Wi, beta, alphaG, epsilon, u, v, psi, wz, Txx, Txy, Tyy
```

This data could feed:

- supervised regression;
- error models;
- surrogates;
- neural operators;
- regime classifiers;
- hybrid physics/statistics models.

## Quality rule

Every port must follow the cycle:

```text
read Fortran -> document the function -> implement Python -> test -> compare -> version
```

No heroic late-night transcriptions without tests. That is like solving
Navier-Stokes on a napkin: romantic, but dangerous.
