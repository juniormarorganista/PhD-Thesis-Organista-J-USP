"""Small analytical manufactured field used for smoke tests.

This is not yet one of the thesis Fortran cases.  It exists only to validate the
Python package structure, grid conventions and derivative/error utilities before
the faithful port of ``case_01.f`` ... ``case_14.f``.
"""

from __future__ import annotations

import numpy as np

from mms_viscoelastic.fields import Field2D
from mms_viscoelastic.grids import Grid2D


def evaluate_toy_solution(grid: Grid2D) -> Field2D:
    """Evaluate a smooth divergence-free toy field on ``grid``."""
    x = grid.X
    y = grid.Y
    psi = np.sin(np.pi * x) * np.sin(np.pi * y)
    u = np.pi * np.sin(np.pi * x) * np.cos(np.pi * y)
    v = -np.pi * np.cos(np.pi * x) * np.sin(np.pi * y)
    wz = 2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    txx = 1.0 + 0.1 * np.sin(np.pi * x)
    txy = 0.1 * np.sin(np.pi * x) * np.sin(np.pi * y)
    tyy = 1.0 + 0.1 * np.cos(np.pi * y)
    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)
