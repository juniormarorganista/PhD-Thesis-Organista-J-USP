"""Reference-data comparison utilities for the Case 12 exact fields.

This module compares the Python manufactured-solution fields
(``evaluate_case_12_exact``) against CSV files generated from an isolated
Fortran reference program that ports ``exact_solutions_case_12(t)``.

Scope
-----
The comparison covers the **t=0 initial-instant** exact fields ``u``, ``v``,
``wz``, ``psi``, ``txx``, ``txy`` and ``tyy``. It does not reproduce the
time-marching solver loop (``source_term_case_12`` and the discretized PDE
update). See ``docs/case-12-solver-comparison.md`` for the full scope
statement.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.validation.source_term_reference import (
    SourceTermReferenceComparison,
    compare_case12_reference_csv,
)

CASE12_FIELD_COMPONENTS: tuple[str, ...] = ("u", "v", "wz", "psi", "txx", "txy", "tyy")


def compare_case12_field_reference(
    path: str | Path,
    *,
    component: str,
    t: float = 0.0,
    parameters: Case12Parameters | None = None,
    atol: float = 1.0e-10,
    rtol: float = 1.0e-10,
) -> SourceTermReferenceComparison:
    """Compare a Python Case 12 exact field against a Fortran-generated CSV.

    Parameters
    ----------
    path:
        Reference CSV with columns ``x,y,u,v,wz,psi,txx,txy,tyy``.
    component:
        Field component to compare. Must be one of ``CASE12_FIELD_COMPONENTS``.
    t:
        Physical time. The reference program defaults to ``t=0``.
    parameters:
        Case-12 parameters. If omitted, defaults match ``Case12Parameters()``,
        which is also the default used by the Fortran reference program.
    """

    if component not in CASE12_FIELD_COMPONENTS:
        raise ValueError(
            f"component must be one of {CASE12_FIELD_COMPONENTS!r}, got {component!r}."
        )

    params = parameters or Case12Parameters()

    def evaluator(grid: Grid2D) -> np.ndarray:
        fields = evaluate_case_12_exact(grid, t=t, parameters=params)
        value = getattr(fields, component)
        if value is None:
            raise ValueError(f"Field component {component!r} was not computed.")
        return value

    return compare_case12_reference_csv(
        path,
        component=component,
        evaluator=evaluator,
        atol=atol,
        rtol=rtol,
    )


__all__ = [
    "CASE12_FIELD_COMPONENTS",
    "compare_case12_field_reference",
]
